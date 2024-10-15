function thr = kthselect(x,tptr,lev,wname)
% 经 Mr.看海 修改后的小波阈值滤波阈值选择函数（原函数为thselect）
% 加入了 VisuShrink 方法，即通用阈值去噪方法
% 该方法的实现主要参考：https://core.ac.uk/download/pdf/217399263.pdf （Denoising techniques - a comparison）
% 如果要使用VisuShrink方法，需要令TPTR =
% 'visushrink'，注意大小写，并指定lev和wname值。此时输出thr为lev*1维度的数组
% 其中lev：  小波分解水平
%     wname：小波名称

%  本代码为淘宝买家专用，不开源，请勿公开分享~


%THSELECT Threshold selection for de-noising.
%   THR = THSELECT(X,TPTR) returns threshold X-adapted value
%   using selection rule defined by string TPTR.
%
%   Available selection rules are:
%   TPTR = 'rigrsure', adaptive threshold selection using
%       principle of Stein's Unbiased Risk Estimate.
%   TPTR = 'heursure', heuristic variant of the first option.
%   TPTR = 'sqtwolog', threshold is sqrt(2*log(length(X))).
%   TPTR = 'minimaxi', minimax thresholding.
%
%   Threshold selection rules are based on the underlying
%   model y = f(t) + e where e is a white noise N(0,1).
%   Dealing with unscaled or nonwhite noise can be handled
%   using rescaling output threshold THR (see SCAL parameter
%   in WDEN).
%
%   See also WDEN.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 20-Dec-2010.
%   Copyright 1995-2020 The MathWorks, Inc.
%narginchk(1,2);
validateattributes(x, {'numeric'}, {'2d'}, 'THSELECT', 'X')
if isStringScalar(tptr)
    tptr = convertStringsToChars(tptr);
end

if isrow(x)
    x = x(:);
end
[n,m] = size(x);

switch tptr
    case 'rigrsure'
        sx = sort(abs(x),1);
        sx2 = sx.^2;
        N1 = repmat((n-2*(1:n))',1,m);
        N2 = repmat((n-1:-1:0)',1,m);
        CS1 = cumsum(sx2,1);
        risks = (N1+CS1+N2.*sx2)./n;
        [~,best] = min(risks,[],1, 'linear');
        % thr will be row vector
        thr = sx(best);
        
    case 'heursure'
        %%
        hthr = (2*log(n)).^0.5;
        eta = sum(abs(x).^2-n,1)./n;
        crit = (log(n)/log(2))^(1.5)/(n.^0.5);
        thr = thselect(x,'rigrsure');
        thr(thr > hthr) = hthr;
        thr(eta < crit) = hthr;
        
    case 'sqtwolog'
        thr = (2*log(n)).^0.5;
        thr = repelem(thr,size(x,2));
        
    case 'minimaxi'
        if n <= 32
            thr = repelem(0, m);
        else
            t = 0.3936 + 0.1829*(log(n)/log(2));
            thr = repelem(t,m);
        end
    case 'visushrink'
        [c,l]=wavedec(x,lev,wname);  
        for i = 1:lev
            d{i} = detcoef(c,l,i); % 提取第i层的高频系数 
            thr(i) = median(abs(d{i}))/0.6745*sqrt(2*log(length(d{i})));
        end
    otherwise
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end


