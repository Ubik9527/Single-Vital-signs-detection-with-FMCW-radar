function y = kwthresh(x,sorh,t,options)
% 经 Mr.看海 修改过后的阈值函数（原函数为wthresh）
% 添加了改进阈值函数类型
% 输入：
% x：阈值函数的输入值，不需修改
% sorh：阈值函数类型
%       sorh ='a1'时，采用改进方法1，改进后的软阈值：
%                     参考论文：孙万麟, 王超. 基于改进的软阈值小波包网络的电力信号消噪[J]. 海军工程大学学报, 2019(4).
%       sorh = 'a2'时，采用改进方法2：
%                     参考论文：刘冲, 马立修, 潘金凤,等. 联合VMD与改进小波阈值的局放信号去噪[J]. 现代电子技术, 2021, 44(21):6.
%       sorh = 'a3'时，采用改进方法3：
%                     参考论文：基于改进小波阈值-CEEMDAN算法的ECG信号去噪研究
%       sorh = 'a4'时，采用改进方法4：
%                     参考论文：基于改进小波阈值函数的语音增强算法研究
%                     采用该方法需要输入两个调节因子，分别是alpha和gamma，取alpha>0，0<gamma<1
%       sorh = 'a5'时，采用改进方法5：
%                     参考论文：《基于VMD与改进小波阈值的地震信号去噪方法研究》
%                     采用该方法需要输入两个调节因子，分别是alpha和beta
% t：阈值
% 输出：
% y：阈值函数输出值
% 测试阈值函数图线可以用下述代码：
% figure('color','w');plot(-5:0.1:5,kwthresh(-5:0.1:5,'a1',1));
% 测试阈值函数图线与软硬阈值图线对比图可用下述代码：
% figure('color','w');plot(-5:0.1:5,kwthresh(-5:0.1:5,'a1',1));hold on
% plot(-5:0.1:5,kwthresh(-5:0.1:5,'s',1))
% plot(-5:0.1:5,kwthresh(-5:0.1:5,'h',1))
% legend('改进阈值方法','软阈值','硬阈值')

%  本代码为淘宝买家专用，不开源，请勿公开分享~

%WTHRESH Perform soft or hard thresholding. 
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's')
%   or hard (if SORH = 'h') T-thresholding  of the input 
%   vector or matrix X. T is the threshold value.
%
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft 
%   thresholding is shrinkage.
%
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard
%   thresholding is cruder.
%
%   See also WDEN, WDENCMP, WPDENCMP.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 24-Jul-2007.
%   Copyright 1995-2020 The MathWorks, Inc.

if isStringScalar(sorh)
    sorh = convertStringsToChars(sorh);
end

switch sorh
  case 's'
    tmp = (abs(x)-t);
    tmp = (tmp+abs(tmp))/2;
    y   = sign(x).*tmp;
 
  case 'h'
    y   = x.*(abs(x)>t);
  case 'a1'
        y = (abs(x)>t).*sign(x).*sqrt(abs(x).^2-t^2);
  case 'a2'
    for i = 1:length(x)
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))-2^(t-abs(x(i))));
        else
            y(i) = 0;
        end
    end
  case 'a3'
    for i = 1:length(x)
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))-2*t/(exp((abs(x(i))-t)/t)+1));
        else
            y(i) = 0;
        end
    end
  case 'a4'  %参考 基于改进小波阈值函数的语音增强算法研究
        alpha = options.a4_alpha; %调节因子
        gamma = options.a4_gamma; %调节因子
        for i = 1:length(x)
            if abs(x(i)) > t
                y(i) = sign(x(i))*(abs(x(i))-(t/(1+alpha).*gamma.^(sqrt(x(i).^2-t.^2))));
            else
                y(i) = sign(x(i))*alpha/(1+alpha)*exp(10*(abs(x(i))-t))*abs(x(i));
            end
        end

  case 'a5' %参考 《基于VMD与改进小波阈值的地震信号去噪方法研究》
    for i = 1:length(x)
        alpha = options.a5_alpha; %调节因子
        beta = options.a5_beta; %调节因子
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))- t/power((abs(x(i)).^beta-(abs(t)).^beta+1),1/beta)*1/exp(power(abs(x(i))-abs(t),1/alpha)));
        else
            y(i) = 0;
        end
    end
    
  otherwise
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end
