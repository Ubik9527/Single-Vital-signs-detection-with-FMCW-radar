%%初始化程序
clear all,clc,close all
t1=clock;
 %% 载入噪声信号数据，数据为.mat格式，并且和程序放置在同一个文件夹下
load('signal.mat');%matrix
YSJ= angle_fft_last;
 %% 数据预处理，数据可能是存储在矩阵或者是EXCEL中的二维数据，衔接为一维的，如果数据是一维数据，此步骤也不会影响数据
[c,l]=size(YSJ);
Y=[];
for i=1:c
    Y=[Y,YSJ(i,:)];
end
[c1,l1]=size(Y);
X=[1:l1];
N1=length(Y);fs =50;
f=(0:N1-1)*(fs/N1); 
 %% 绘制噪声信号图像
figure;
plot(X,Y);
xlabel('横坐标');
ylabel('纵坐标');
title('原始信号');
         
 %% 硬阈值处理
lev=4;
xd=wden(Y,'minimaxi','h','one',lev,'db4');%硬阈值去噪处理后的信号序列
figure
plot(X,xd)
xlabel('横坐标');
ylabel('纵坐标');
title('硬阈值去噪处理')
set(gcf,'Color',[1 1 1])      

 %% 软阈值处理
lev=4;
xs=wden(Y,'heursure','h','one',lev,'db4');%软阈值去噪处理后的信号序列
figure
plot(X,xs)
xlabel('横坐标');
ylabel('纵坐标');
title('软阈值去噪处理')
set(gcf,'Color',[1 1 1])

%% 固定阈值后的去噪处理
lev=4;
xz=wden(Y,'sqtwolog','h','one',lev,'db4');%固定阈值去噪处理后的信号序列
figure
plot(X,xz);
xlabel('横坐标');
ylabel('纵坐标');
title('固定阈值后的去噪处理')
set(gcf,'Color',[1 1 1])

%% 频谱比较
f_xs = abs(fft(xs));            
f_Y = abs(fft(Y));   
f_xd = abs(fft(xd));   
f_xz = abs(fft(xz));            
figure;
plot(f(10:130),f_Y(10:130));hold on;
plot(f(10:130),f_xd(10:130));hold on
plot(f(10:130),f_xs(10:130));hold on
plot(f(10:130),f_xz(10:130));
xlabel('频率(Hz)');ylabel('幅度');title('去噪处理前后信号FFT');
legend('去噪前信号','硬阈值去噪','软阈值去噪','固定阈值后的去噪')
%% 计算信噪比SNR
Psig=sum(Y*Y')/l1;
Pnoi1=sum((Y-xd)*(Y-xd)')/l1;
Pnoi2=sum((Y-xs)*(Y-xs)')/l1;
Pnoi3=sum((Y-xz)*(Y-xz)')/l1;
SNR1=10*log10(Psig/Pnoi1);
SNR2=10*log10(Psig/Pnoi2);
SNR3=10*log10(Psig/Pnoi3);
%% 计算均方根误差RMSE
RMSE1=sqrt(Pnoi1);
RMSE2=sqrt(Pnoi2);
RMSE3=sqrt(Pnoi3);
%% 输出结果
disp('-------------三种阈值设定方式的降噪处理结果---------------'); 
disp(['硬阈值去噪处理的SNR=',num2str(SNR1),'，RMSE=',num2str(RMSE1)]);
disp(['软阈值去噪处理的SNR=',num2str(SNR2),'，RMSE=',num2str(RMSE2)]);
disp(['固定阈值后的去噪处理SNR=',num2str(SNR3),'，RMSE=',num2str(RMSE3)]);
t2=clock;
tim=etime(t2,t1);
disp(['------------------运行耗时',num2str(tim),'秒-------------------'])
