% %% 
% clc;
% clear all;
% close all;
% % 加载信号
% load signal.mat;
% s = angle_fft_last;
% plot(s);title('原始信号');
% [thr,sorh,deepapp,crit]=ddencmp('den','wp',s);
% [x,wpt,perf0,perf1]=wpdencmp(s,sorh,4,'db4',crit,thr,deepapp);
% figure(2);
% plot(x);title('去噪信号（3层）');
% 原始信号参数
fs = 1000;                % 采样频率
t = 0:1/fs:1;             % 时间向量
f = 10;                   % 信号频率
x = sin(2*pi*f*t);        % 原始信号
 
% 添加高斯噪声
snr = 10;                 % 信噪比（以dB为单位）
noise = randn(size(x));   % 生成高斯噪声
noise = noise / sqrt(sum(noise.^2)) * sqrt(sum(x.^2)) / (10^(snr/20)); % 根据信噪比调整噪声的幅度
noisySignal = x + noise;  % 带噪声的信号
 
% 小波去噪参数设置
wavelet = 'morlet';       % 小波基函数
level = 5;            % 小波分解的层数
threshold = 'soft';   % 阈值选择方法为 'soft'
 
% 小波去噪
denoisedSignal = wdenoise(noisySignal, level, ...
    'Wavelet', wavelet, 'DenoisingMethod', 'UniversalThreshold', 'ThresholdRule', threshold);
 
% 绘制加噪信号和去噪信号
t = 1:length(noisySignal);
figure;
subplot(2,1,1);
plot(t, noisySignal);
title('加噪信号');
xlabel('时间');
ylabel('幅值');
subplot(2,1,2);
plot(t, denoisedSignal);
title('去噪信号');
xlabel('时间');
ylabel('幅值');
 