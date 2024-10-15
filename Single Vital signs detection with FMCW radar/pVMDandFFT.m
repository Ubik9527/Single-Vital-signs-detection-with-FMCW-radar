function [imf,CenFs] = pVMDandFFT(y,FsOrT, alpha, K, tol)
% 画信号VMD分解与各IMF分量频谱对照图
% 输入：
% y：待分解的数据(一维)
% FsOrT： 采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量。如果未知采样频率，可设置为1
% alpha   - 惩罚因子
% K       - 指定分解模态数
% tol     - 收敛容差，是优化的停止准则之一，可以取 1e-6~5e-6

% 输出：
% imf：内涵模态分量，统一为n*m格式，其中n为模态数，m为数据点数。例如 imf(1,:)即IMF1，imf(end,:)即为残差
% 注意，为了与其他“类EMD”方法分解出来的imf分量保持一致，本程序内将imf排序进行了翻转
% 即保证imf排列是从高频向低频排列
% CenFs：即CentralFrequencies，各imf分量的中心频率

% 例1：（FsOrT为采样频率）
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pVMDandFFT(y,fs,2000,2,1e-7);
% 例2：（FsOrT为时间向量，需要注意此时FsOrT的长度要与y相同）
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pVMDandFFT(y,t,2000,2,1e-7);

% 注意：在使用该代码之前，请务必安装工具箱：http://www.khscience.cn/docs/index.php/2020/04/09/1/

%  Copyright (c) 2021 Mr.括号 All rights reserved.
%  本代码为淘宝买家专用，不开源，请勿公开分享~
%%
if length(FsOrT) == 1
    t = 1/FsOrT:1/FsOrT:length(y)/FsOrT;
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
[imf,CenFs] = kVMD(y, FsOrT,alpha, K, tol);
figure('Name','VMD分解与各IMF分量频谱对照图','Color','white');
subplot(size(imf,1)+1,2,1);
plot(t,y,'k');grid on;
ylabel('原始数据');
title('VMD分解');
set(gca,'XTick',[]);
subplot(size(imf,1)+1,2,2);
pFFT(y,Fs);grid on;
title('对应频谱');
set(gca,'XTick',[]);
for i = 2:size(imf,1)+1
    subplot(size(imf,1)+1,2,i*2-1);
    plot(t,imf(i-1,:),'k');
    ylabel(['IMF',num2str(i-1)]);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        ylabel(['res']);
        xlabel('time/s');
    end
    grid on;
    subplot(size(imf,1)+1,2,i*2);
    pFFT(imf(i-1,:),Fs);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        xlabel('frequency/Hz');
    end
    grid on;
end
end
