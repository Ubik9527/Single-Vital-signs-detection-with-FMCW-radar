%% 通过VMD方法得到心率与呼吸率
% 输入参数：
% x：输入生命体征信号
% fs：信号频率
% 输出参数：
% breathcount：呼吸率
% heartcount：心率
function [breathcount,heartcount] = result_VMD(x,fs)
FsOrT = fs;
alpha = 2000;
K = 3;
tol = 1e-6;
[imf,CenFs] = pVMDandFFT(x,FsOrT, alpha, K, tol) ;

%% VMD得到心率与呼吸率
N = size(imf,2);f=(0:N-1)*(fs/N);    
for i =1:size(imf,1)
    FFT_imf = abs(fft(imf(i,:)));   
    [pk(i),loc(i)] = max(FFT_imf);
end
heart_f = [];
if ((loc(2)-1)*fs/N >0.85 && (loc(2)-1)*fs/N <1.8)
    heart_f(end+1) = loc(2);
end
if ((loc(1)-1)*fs/N >0.85 && (loc(1)-1)*fs/N <1.8)
    heart_f(end+1) = loc(1);
end
if (length(heart_f)>1)
    if(pk(1)>pk(2))
        heartcount = (loc(1)-1)*fs/N*60;
    else 
        heartcount = (loc(2)-1)*fs/N*60;
    end
else
    heartcount = (heart_f-1)*fs/N*60;
end
if((loc(3)-1)*fs/N >0.15 && (loc(3)-1)*fs/N <0.4)
    breathcount = (loc(3)-1)*fs/N*60;
elseif ((loc(2)-1)*fs/N >0.15 && (loc(2)-1)*fs/N <0.4) 
    breathcount = (loc(2)-1)*fs/N*60;
else
    breathcount = 0;
end

    % if((loc(4)-1)*fs/N >0.15 && (loc(4)-1)*fs/N < 1.8)
    %     breathcount = (loc(4)-1)*fs/N*60;
    % else
    %     breathcount = (loc(3)-1)*fs/N*60;
    % end
    % if((loc(3)-1)*fs/N >0.8 && (loc(3)-1)*fs/N <2.0)
    %     heartcount = (loc(3)-1)*fs/N*60;
    % elseif ((loc(2)-1)*fs/N >0.8 && (loc(2)-1)*fs/N <1.8)
    %         heartcount = (loc(2)-1)*fs/N*60;
    % else heartcount = (loc(1)-1)*fs/N*60;
% end
end