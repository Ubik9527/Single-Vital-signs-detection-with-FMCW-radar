function pFFT(y,Fs)
%  该函数有两个输入，即信号值与采样频率
%  该函数只画一张图，内置plot，但未使用figure
%  Copyright (c) 2020 Mr.括号 All rights reserved.
%  本代码为淘宝买家专用，不开源，请勿公开分享~
%% 
t_s = 1/Fs; %采样周期
t_start = 0;
t_end = (length(y)-1)*t_s;
t = 0 : t_s : t_end;
Druation = t_end -t_start;  %计算采样时间
Sampling_points = Druation/t_s +1;  %采样点数
nfft = Sampling_points;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%频谱%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_f = fft(y,nfft); %傅里叶变换
f_s = 1/t_s; %采样频率
f_x = 0:f_s/(Sampling_points -1):f_s;  %注意这里和横坐标频率对应上了，频率分辨率就是f_s/(Sampling_points -1)
plot(f_x(1:round(length(f_x)/2)),2/Sampling_points*abs(y_f(1:round(length(f_x)/2))),'k');
