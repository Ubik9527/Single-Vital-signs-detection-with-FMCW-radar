%% 硬件平台 ：AWR1843+DCA1000EVM
%% 功能：单人呼吸心跳原始数据采集与生命体征信号处理与提取
%% 采集环境：办公桌前，胸部正对采集设备0.5m距离
%% 算法流程：预处理+MTI +反正切提取相位+相位解缠+相位差分+滑动平均去噪+带通滤波器
%% ========================================================================
clc;
clear;
close all;

%% 雷达参数设置
numADCSamples = 256; % number of ADC samples per chirp采样点数200
numADCBits = 16;     % number of ADC bits per sample
numTX = 1;                %发射天线数
numRX = 4;           % number of receivers：接收天线数
numLanes = 2;        % do not change. number of lanes is always 2
isReal = 0;          % set to 1 if real only data, 0 if complex data：1为实采样，0为复采样
chirpLoop = 2;
numframes = 3000;
%1T4R（启用一条发射天线TX0）256个采样点，共3000帧，帧周期20ms，共60s
Fs=6e6;             %ADC采样率 ：6Msps            
c=3*1e8;            %光速
ts=numADCSamples/Fs;    %ADC采样时间
slope=70e12;                 %调频斜率 ：~70.006MHz/us
B_valid =ts*slope;                %有效带宽：3GHz
detaR=c/(2*B_valid);          %距离分辨率：5cm
startFreq = 77e9;              %起始频率 ：77GHz
lambda=c/startFreq;              % 雷达信号波长：  λ=c/f0 = 5mm
virtualAntenna=numRX*numTX;
parameter  = generateParameter();
fs = 50;

%% 读取传感器的参考数据
addpath(genpath('C:\Users\Lee\Desktop\radar data'));
%% 读取呼吸传感器数据
Filename1 = '20240131144003.txt';
fid = fopen(Filename1,'r');
inputdata_heart=  fscanf(fid,'%x');
fclose(fid);
filesize = size(inputdata_heart, 1);
%筛选最后两位
for i=1:(filesize/7)
    rawdata_breath(i,:) = inputdata_heart((i-1)*7+1:i*7);
end
%合成16位
rawdata_breath = bitshift(rawdata_breath(:,6),8)+rawdata_breath(:,7);
rawdata_breath1 = rawdata_breath(1:3000);
%去除直流分量l
rawdata_breath = rawdata_breath1 - mean(rawdata_breath1);

%% 心率带数据读取
% % 读取图像文件
% image = imread('tfy_static_7.jpg'); % 这里的'your_image.jpg'应该被替换为你要读取的图像文件路径及名称
%  
% % 显示图像
% figure;
% imshow(image);
% title('心率数据'); % 设置标题（可选）
%% 读取心跳传感器数据
Filename2 = '20240131144004.txt';
%读取原始数据
fid = fopen(Filename2,'r');
inputdata_heart=  fscanf(fid,'%x');
fclose(fid);
filesize = size(inputdata_heart, 1);
%筛选最后两位
for i=1:(filesize/7)
    rawdata_heart(i,:) = inputdata_heart((i-1)*7+1:i*7);
end
%合成16位
rawdata_heart = bitshift(rawdata_heart(:,6),8)+rawdata_heart(:,7);
for i=1:fix((size(rawdata_heart)/4))
     data_average(i)=(rawdata_heart(4*i-3)+rawdata_heart(4*i-2)+rawdata_heart(4*i-1)+rawdata_heart(4*i))/4;
end
rawdata_heart1 = data_average(51:3050);
%去除直流分量
rawdata_heart = rawdata_heart1 - mean(rawdata_heart1);

%% 传感器测得呼吸心跳频率
% N = length(rawdata_breath);
% heart_raw = abs(fft(rawdata_breath));
% [index_rawbreath,~] = findmax(heart_raw,5,30);
% rr_raw = (index_rawbreath-1)*fs/N*60;
% heart_raw = abs(fft(rawdata_heart));
% [index_rawbreath,~] = findmax(heart_raw,50,120);
% hr_raw = (index_rawbreath-1)*fs/N*60;
% disp(['传感器测得结果： 呼吸：',num2str(rr_raw),'  心跳：',num2str(hr_raw)])
%% 读取Bin文件（数据预处理）
Filename = 'tfy_static_5_Raw_0.bin';       
fid = fopen(Filename,'r');
adcDataRow = fread(fid, 'int16');
if numADCBits ~= 16 
    l_max = 2^(numADCBits-1)-1;
    adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - 2^numADCBits;
    %16bit的ADC，最高表示数据为[-32768，32767]，即[-2^15,2^15]，
    % 剩下的1位的符号位，程序中数值高于32767的都要减去65536。
end
fclose(fid);
process_num = numframes;                     %只对process_num帧的数据做处理
%  fileSize1 = size(adcDataRow, 1);
fileSize=process_num*numADCSamples*numTX*numRX*2;   
PRTnum = fix(fileSize/(numADCSamples*numRX));
fileSize = PRTnum * numADCSamples*numRX;
adcData = adcDataRow(1:fileSize);
% real data reshape, filesize = numADCSamples*numChirps
if isReal
    numChirps = fileSize/numADCSamples/numRX;
    LVDS = zeros(1, fileSize);
    %create column for each chirp
    LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
else%复采样
    numChirps = fileSize/2/numADCSamples/numRX;     %含有实部虚部，除以2
    %共2048个chirps（1024帧*2个chirp）
    LVDS = zeros(1, fileSize/2);
    %combine real and imaginary part into complex data将实部虚部结合成复数
    %read in file: 2I is followed by 2Q     adcData数据组成:两个实部，接着是两个虚部
    counter = 1;
    for i=1:4:fileSize-1        %1T4R
        LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2);       %复数形式
        LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
        counter = counter + 2;
    end
    % create column for each chirp：每一列为chirp
    LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
    %each row is data from one chirp：每一行为chirp
    LVDS = LVDS.';
end

%% 重组数据（4条接收天线的复数数据）
adcData = zeros(numRX,numChirps*numADCSamples);
for row = 1:numRX
    for i = 1: numChirps
        adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
    end
end

%% 3DFFT
rangeRes     = detaR; %距离分辨率 有效带宽
rangeIndex   = (0:numADCSamples-1) * rangeRes;
speedRes     = lambda / (2 * numframes* ts);
dopplerIndex = (-numframes/2:1:numframes/2 - 1) * speedRes;
rawData      = reshape(adcData,virtualAntenna,numADCSamples, numChirps);
channelNum    = size(rawData,1);
rangebinNum   = size(rawData,2);
frameNum      = size(rawData,3);
N=rangebinNum;M=frameNum;Q=1024;
adcData2=adcData.';
rawdata_radar=reshape(adcData2,rangebinNum,frameNum,channelNum);
% 距离FFT
range_win = hamming(rangebinNum);   %加海明窗
doppler_win = hamming(frameNum);
range_profile = [];
for k=1:channelNum
   for m=1:frameNum
      temp=rawdata_radar(:,m,k).*range_win;    %加窗函数
      temp_fft=fft(temp,N);    %对每个chirp做N点FFT
      range_profile(:,m,k)=temp_fft;
   end
end

% % 多普勒FFT
% speed_profile = [];
% for k=1:channelNum
%     for n=1:N
%       temp=range_profile(n,:,k).*(doppler_win)';    
%       temp_fft=fftshift(fft(temp,M));    %对rangeFFT结果进行M点FFT
%       speed_profile(n,:,k)=temp_fft;  
%     end
% end

%% 静态滤波
% 选中一个通道的数据
data_input = range_profile(:,:,1);
data_out = static_filtering(data_input,2,rangeRes);
% 静态滤波方法选择：
%   0：不滤波
%   1：脉冲对消法
%   2：平均相消法 *
%   3：相位判断法


%% 生命体征信号提取
t = 5;
[vital_sign,range_index] = ranngebin_select(data_out,fs,t,rangeRes,4);
% n：距离门选择方法：
%     1：固定距离门
%     2：最大距离门
%     3：根据一段时间内的距离门单元众数选择
%     4：根据一段时间内的非相干累积值选择


%% 进行相位解缠
%或称为相位解卷绕，由于相位值在 [ − π ,π ] 之间，而我们需要相位展开以获取实际的位移曲线，
%因此每当连续值之间的相位差大于或者小于±π时，通过从相位中减去2π来获得相位展开。
%unwrap函数：
phi=unwrap(vital_sign); 
figure;
plot(1/fs:1/fs:M/fs,phi);
title('相位解缠绕后信号');
xlabel('时间/s')
ylabel('相位/rad')

%% 不进行运动检测直接计算呼吸率心率
[breath_count,heart_count,breath_data,heart_data,breath_spectrum,heart_spectrum] = spectral_method(phi,fs);
breath_directdata =breath_data;
breath_direct = breath_spectrum;
heart_directdata = heart_data;
heart_direct = heart_spectrum;
disp(['分段距离门结果 呼吸：',num2str(breath_count),'  心跳：',num2str(heart_count)])
%信噪比 
% SNR_heart = 10*log10(heart(heart_index)^2/(sum(breath.^2)-heart(heart_index)^2));
% disp(['的那段信噪比结果 呼吸：',num2str(SNR_breath),'  心跳：',num2str(SNR_heart)])

%% SARIMA运动修复
phi_recover = motionrecovery(phi,0.2);

%% phase difference 相位差分
phi_recover = phi;
for i = 2:length(phi_recover)
    phi_recoverdiff(i) = phi_recover(i) - phi_recover(i-1);
end 
% 脉冲噪声去除
for i=2:length(phi_recoverdiff)
    if(abs(phi_recoverdiff(i))>5*mean(abs(phi_recoverdiff)))
        phi_recoverdiff(i)=phi_recoverdiff(i-1);
    end
end
FFT_raw = abs(fft(phi_recoverdiff));

%% 去除呼吸谐波
%先找到符合呼吸频率区间幅度最大的频率
[fft_index,~ ] = findmax(FFT_raw,5,25);
% 二次谐波去除
N1 = length(phi_recoverdiff);
f=(0:N1-1)*(fs/N1);  
breath_harmonic2=fs/N1*(fft_index-1)*2;
phi_recoverdiff2 = filter(bandstop(breath_harmonic2),phi_recoverdiff);
% 三次谐波去除
breath_harmonic3=fs/N1*(fft_index-1)*3;
phi_bandstop = filter(bandstop(breath_harmonic3),phi_recoverdiff2);
FFT_bandstop = abs(fft(phi_bandstop));  
figure;
plot(f(1:N1/8),FFT_raw(1:N1/8));hold on;
plot(f(1:N1/8),FFT_bandstop(1:N1/8))
xlabel('频率(Hz)');ylabel('幅度');title('带通滤波前后'); legend('原信号','带阻滤波后');

%% 小波阈值去噪
% xb=wden(phi_bandstop,'heursure','s','sln',4,'sym10');
%% VMD+小波阈值
xb = filEMDsWaveletTh(phi_bandstop);
xb_f = abs(fft(xb));
%傅里叶变换结果对称
figure;
plot(f(1:N1/8),FFT_bandstop(1:N1/8));hold on;
plot(f(1:N1/8),xb_f(1:N1/8)); 
xlabel('频率（Hz）');ylabel('幅度');title('相位信号FFT');legend('去噪前','小波去噪后');

%%  IIR4阶带通滤波 Bandpass Filter 0.1-0.5hz，输出呼吸信号
breath_data = filter(IIR_breath,xb); 
breath = abs(fft(breath_data));       
figure;
X=(1:size(breath_data,2))/fs;
plot(X,breath_data/abs(max(breath_data)));hold on;
plot(X,breath_directdata/abs(max(breath_directdata)));
xlabel('时间(s)');ylabel('幅度(归一化)');title('呼吸时域波形');legend('运动恢复加小波滤波后波形','原始波形'); 
figure;
plot(f(1:130),breath(1:130));hold on;  %取前一部分放大观察
plot(f(1:130),breath_direct(1:130));
xlabel('频率(Hz)');ylabel('幅度');title('呼吸信号FFT');legend('运动恢复加小波滤波后频谱','原始频谱');
[breath_index,~] = findmax(breath,5,25);
breath_count_mix =60*(breath_index-1)*(fs/N1);
SNR_breath = 10*log10(breath(breath_index)^2/(sum(breath.^2)-breath(breath_index)^2));
%% IIR8阶带通滤波 Bandpass Filter 0.8-2hz，输出心跳的数据
heart_data = filter(IIR_heart,xb); 
heart = abs(fft(heart_data));
figure;
X=(1:size(heart_data,2))/fs;
plot(X,heart_data/abs(max(heart_data)));hold on;
plot(X,heart_directdata/abs(max(heart_directdata)));
xlabel('时间(s)');ylabel('幅度(归一化)');title('心跳时域波形');legend('运动恢复加小波滤波后波形','原始波形');
figure;
plot(f(1:130),heart(1:130));hold on;
plot(f(1:130),heart_direct(1:130));
xlabel('频率(Hz)');ylabel('幅度');title('心跳信号FFT');legend('运动恢复加小波滤波后频谱','原始频谱');
[heart_index,max_value] = findmax(heart,51,110);
heart_count_mix =60*(heart_index-1)*(fs/N1);
SNR_heart = 10*log10(heart(heart_index)^2/(sum(heart.^2)-heart(heart_index)^2));
%% 频谱法测量结果
disp(['小波滤波恢复结果 呼吸：',num2str(breath_count_mix),'  心跳：',num2str(heart_count_mix)])
%% 短时傅里叶测量结果
wlen = 1024;nfft = 1024;hop = wlen/4;fs = 50;x=heart_data;
win = blackman(wlen, 'periodic');
[S, f, t] = mystft(x, win, hop, nfft, fs);
for i=1:length(t)
    [~,fmaxindex(i)] = max(S((wlen/2+1):wlen,i));
    fmax(i) = f(wlen/2+fmaxindex(i))*60;
end
majority_value = mode(fmax);
average_value = mean(fmax);
median_value = median(fmax);
disp(['短时傅里叶心率结果',' 平均数：',num2str(average_value),' 众数：',num2str(majority_value),' 中位数：',num2str(median_value)])

%% 峰值检测法输出结果
[pks,~] = findpeaks(breath_data,'MinPeakDistance',125);
breath_peakcount = length(pks);
[pks,~] = findpeaks(heart_data,'MinPeakDistance',20);
heart_peakcount = length(pks);
disp(['峰值检测法结果 呼吸：',num2str(breath_peakcount),'  心跳：',num2str(heart_peakcount)])

%% VMD测量结果
[breath_vmdcount,heart_vmdcount] = result_VMD(xb,fs);
disp(['VMD方法结果 呼吸:',num2str(breath_vmdcount),'  心率： ',num2str(heart_vmdcount)]);

%% 提取多个峰值
% [pks_breath,locs_breath] = findpeaks(breath(5:30),'SortStr','descend');
% [pks_heart,locs_heart] = findpeaks(heart(50:120),'SortStr','descend');
% for i=1:min(length(locs_breath),5)
%     disp(['呼吸第',num2str(i),'可能：',num2str((locs_breath(i)-1)*fs/N1*60),'   归一化幅值：',num2str(pks_breath(i)/pks_breath(1))])
% end
% for i=1:min(length(locs_heart),10)
%     disp(['心跳第',num2str(i),'可能：',num2str((locs_heart(i)-1)*fs/N1*60),'   归一化幅值：',num2str(pks_heart(i)/pks_heart(1))])
% end
% %信噪比 
% SNR_heart = 10*log10(heart(heart_index)^2/(sum(breath.^2)-heart(heart_index)^2));

%输出运动失真恢复后信号
% disp(['运动失真恢复后信噪比结果 呼吸：',num2str(SNR_breath),'  心跳：',num2str(SNR_heart)])


