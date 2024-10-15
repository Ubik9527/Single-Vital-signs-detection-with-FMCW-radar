function [imf,CenFs] = pVMD(y,FsOrT, alpha, K, tol)
% 画信号VMD分解图
% 输入：
% data：待分解的数据(一维)
% FsOrT：采样频率或采样时间向量，如果为采样频率，该变量输入单个值；如果为时间向量，该变量为与y相同长度的一维向量。如果未知采样频率，可设置为1
% alpha   - 惩罚因子
% K       - 指定分解模态数
% tol     - 收敛容差，是优化的停止准则之一，可以取 1e-6~5e-6

% 输出：
% imf：内涵模态分量，统一为n*m格式，其中n为模态数，m为数据点数。例如 imf(1,:)即IMF1，imf(end,:)即为残差
% 注意，为了与其他“类EMD”方法分解出来的imf分量保持一致，本程序内将imf排序进行了翻转
% 即保证imf排列是从高频向低频排列
% CenFs：即CentralFrequencies，各imf分量的中心频率

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
%% 1.绘制二维图
figure('Name','VMD分解与各IMF分量对照图','Color','white');
subplot(size(imf,1)+1,1,1);
plot(t,y,'k');grid on;
ylabel('原始数据');
title('VMD分解');
set(gca,'XTick',[]);

for i = 2:size(imf,1)+1
    subplot(size(imf,1)+1,1,i);
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

end

%% 2.绘制三维图
figure('Name','模态分解与各IMF分量对照图（三维）','Color','white');
imfall = [y(:)'; imf];
[m, n] = size(imfall);
[X, Y] = meshgrid(1 : m, 1 : n );    % 建立xy平面坐标网格，X、Y都是n行m列的矩阵，X的每一行都是(1:m)，Y的每一列都是(n:-1:1)
Z = imfall(1 : m, :);                    % 选择要显示的时域信号数据，Z是m行n列的矩阵

plot3(X, Y, Z); grid on;

set(gca,'Ydir','reverse'); % 反转y轴
xticks(1:m)
xticklabels(['原始信号', strcat('IMF', string(1:m-2)), 'res']); % 设置x轴的刻度标签

ylabel('数据序列');
% ylabel('时间(s)');
% zlabel('幅值(W)');
zlabel('幅值');
set(gcf,'Color', 'white');

title('VMD分解-三维展开图');

view(-35, 30);
end
