%% 频谱法得到呼吸率和心率结果
% 输入参数：
% angle_data：解缠后的生命体征信号
% fs：采样频率
% 输出参数：
% Breathcount：呼吸率
% Heartcount:心率
% breath_data：呼吸时域波形
% heart_data：心跳时域波形
% breath_spectrum：呼吸频谱
% heart_spectrum：时域频谱

function [breathcount,heartcount,breath_data,heart_data,breath_spectrum,heart_spectrum] = spectral_method(angle_data,fs)
% 相位差分
phi = angle_data;
for i = 2:length(angle_data)
    phi(i) = angle_data(i) - angle_data(i-1);
end 
% 脉冲噪声去除
for i=2:length(phi)
    if(abs(phi(i))>5*mean(abs(phi)))
        phi(i)=phi(i-1);
    end
end
% 滑动平均滤波
phi=smoothdata(phi,'movmean',5);
%对相位信号作FFT 
N=length(phi);
FFT = abs(fft(phi));            %--FFT   取模，幅度
f=(0:N-1)*(fs/N);             %其中每点的频率
% 去除呼吸谐波
fft_max = 0; 
for i = 5:25
    if (FFT(i) > fft_max)    
        fft_max = FFT(i);
        if(fft_max<1e-2)          %幅度置信 判断是否是存在人的心跳
            fft_index=3001;%不存在
        else
            fft_index=i;
        end
    end
end

% 带阻滤波器
% 二次谐波去除
breath_harmonic2=fs/N*(fft_index-1)*2;
phi_bandstop2 = filter(bandstop(breath_harmonic2),phi);
% 三次谐波去除
breath_harmonic3=fs/N*(fft_index-1)*3;
phi_bandstop3 = filter(bandstop(breath_harmonic3),phi_bandstop2);
phi=phi_bandstop3;
% FFT_bandstop = abs(fft(phi_bandstop));  
% figure;
% plot(f(1:N1/8),FFT2(1:N1/8)) %取前一部分放大观察
% xlabel('频率（f/Hz）');
% ylabel('幅度');
% title('带阻滤波后原始相位信号FFT ');
%  IIR带通滤波呼吸信号 
breath_data = filter(IIR_breath,phi); 
% 谱估计
N=length(breath_data);
breath_spectrum = abs(fft(breath_data));            
breath_max = 0;
%谱峰最大值搜索
for i =1:length(breath_spectrum)
    if ((i-1)*fs/N > 0.15) && ((i-1)*fs/N <0.4)
        if(breath_spectrum(i)>breath_max)
            breath_max = breath_spectrum(i);
            breath_index=i;
        end
    end
end
%呼吸频率解算
breathcount =60*(breath_index-1)*(fs/N);
%计算信噪比
SNR_breath = 10*log10(breath_spectrum(breath_index)^2/(sum(breath_spectrum.^2)-breath_spectrum(breath_index)^2));
% IIR带通滤波心跳信号 
heart_data = filter(IIR_heart,phi); 
N=length(heart_data);
heart_spectrum=abs(fft(heart_data));
heart_max = 0;
%谱峰最大值搜索
for i =1:length(breath_spectrum)
    if ((i-1)*fs/N > 0.82) && ((i-1)*fs/N <2.0)
        if(heart_spectrum(i)>heart_max)
            heart_max = heart_spectrum(i);
            heart_index=i;
        end
    end
end
%呼吸频率解算
heartcount =60*(heart_index-1)*(fs/N);
SNR_heart = 10*log10(heart_spectrum(heart_index)^2/(sum(heart_spectrum.^2)-heart_spectrum(heart_index)^2));
end