function reSig = filEMDsWaveletTh(sig)
% ICEEMDAN-样本熵-改进的小波阈值的封装函数，该函数中的ICEEMDAN方法、样本熵、小波阈值的相关参数可以进行修改替换
% 输入参数:
%   sig - 待滤波信号
%
% 输出参数:
%   reSig - 滤波后的信号

%% 1.ICEEMDAN分解
% Nstd = 0.2; %Nstd为附加噪声标准差与Y标准差之比
% NE = 100;   %NE为对信号的平均次数
% MaxIter = 1000;% MaxIter 最大迭代次数
% imf = pICEEMDAN(sig,1,Nstd,NE,MaxIter);
FsOrT = 50;
alpha = 2000;
K = 6;
tol = 1e-6;
imf = pVMDandFFT(sig,FsOrT, alpha, K, tol);
%% 2.熵特征提取，如果要更换其他熵特征，可以参考testGenFeaEn.m文件的演示
featureNamesCell = {'SpEn'}; %要进行特征提取的特征名称，可以根据实际需要进行删减，留下的特征注意拼写正确

% 设置样本熵参数，Spdim为样本熵的模式维度，Spr为样本熵的阈值，如果不提取样本熵特征可以删除以下两行
option.Spdim  = 2;
option.Spr   = 0.15;

fea = genFeatureEn(imf,featureNamesCell,option);  %调用genFeature函数，完成特征提取，算出的特征值会保存在fea变量里，
                                             %fea变量的长度和featureNamesCell中指定的特征量一致，且顺序一一对应
                                             %程序运行完成后，在MATLAB的工作区，双击fea变量，可以查看求得的具体数值
% disp(['样本熵值为:',num2str(fea')])
% % 绘制熵数值图，用于分析各IMF分量的复杂性
% figure('Color','w');
% bar(fea);  % 生成条形图显示每个IMF分量的多尺度排列熵
% hold on;
% plot(fea, 'r--o');  % 在条形图上覆盖绘制折线图，以便直观比较
% title(featureNamesCell{1});
% xlabel('IMF');
% ylabel('熵值');

% %% 3.根据熵值筛选IMF分量
% th = 0.1; %阈值
% selectedIMF = imf(fea>th,:); % 筛选样本熵大于0.1的IMF分量
% selectedIdx = find(fea>th); % 找出满足条件的IMF分量编号
% selectedIMF = imf(1:5,:);
% selectedIdx = [1,2,3,4,5].';
% % fprintf('满足条件的IMF分量编号为:');
% % disp(selectedIdx');

%% 3. 根据IMF与原信号的相关性筛选IMF分量
thre = 0.5; %阈值
for i = 1:K
    R = corrcoef(sig,imf(i,:));
    r(i) = R(1,2);
end
fprintf('IMF分量的相关系数:');
disp(r);
selectedIMF = imf(r<thre,:); % 筛选样本熵大于0.1的IMF分量
selectedIdx = find(r<thre); % 找出满足条件的IMF分量编号
fprintf('需要滤波的IMF分量为:');
disp(selectedIdx);
%% 3.根据频率筛选IMF分量
% selectedIdx = [];
% selectedIMF = [];
% N = size(imf,2);fs = 50;   
% for i = 1:size(imf,1)
%     [~,index] =  max(abs(fft(imf(i,:))))
%     main_f = (index-1)*fs/N;
%     if(main_f > 0.1 && main_f <2.0)
%         selectedIdx(end+1) = main_f;
%         selectedIMF(i,:) = imf(i,:);
% end
% fprintf('满足条件的IMF分量编号为:');
% disp(selectedIdx);
%% 4.对选出的IMF分量进行小波阈值滤波
wname = 'sym10'; % 小波名称
SORH = 'a3'; % 阈值函数,采用改进的软阈值a4
lev = 4; % 分解层数
tptr = 'heursure'; % 阈值选择规则
options.a4_alpha = 2; % a4方法的alpha参数
options.a4_gamma = 0.8; % a4方法的gamma参数

denoisedIMFs = zeros(size(selectedIMF)); % 初始化

N_selected = size(selectedIMF,1);
figure('Color','w');
N_selected = size(selectedIMF,1);
for i=1:N_selected
    denoisedIMFs(i,:) = filterWaveletTh(selectedIMF(i,:),wname,SORH,lev,tptr,options);
    % 绘制滤波前后对比图
    subplot(N_selected,2,2*i-1);
    plot(selectedIMF(i,:),'k');
    title(['IMF ',num2str(selectedIdx(i)),' (滤波前)']);
    ylim_curr = ylim; % 记录当前子图的y轴范围
    
    subplot(N_selected,2,2*i);
    plot(denoisedIMFs(i,:),'k');
    title(['IMF ',num2str(selectedIdx(i)),' (滤波后)']);
    ylim(ylim_curr); % 将y轴范围设置为与左侧子图相同
end
%% 5.重构滤波后的信号
reSig = sum(denoisedIMFs,1);

% 将未筛选的IMF分量直接相加到结果中
unselectedIMF = imf(~ismember(1:size(imf,1),selectedIdx),:);
reSig = reSig + sum(unselectedIMF,1);

end