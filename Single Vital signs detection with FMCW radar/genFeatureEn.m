function fea = genFeatureEn(data,featureNamesCell,options)
% 熵相关算法的信号特征提取函数
% 输入：
% data：待特征提取的时域信号，可以是二维数据，维度为m*n，其中m为数据组数，n为每组数据的长度。即每行数据为一组。行列方向不可出错
% options：其他设置，使用结构体的方式导入。目前可设置变量包括：
%   -svdpEn：即奇异值的窗口长度。
%   -Apdim：近似熵参数，Apdim为近似熵的模式维度
%   -Apr：近似熵参数，Apr为近似熵的阈值
%   -Spdim：样本熵参数，Spdim为样本熵的模式维度
%   -Spr：Spr为样本熵的阈值
%   -Fuzdim：模糊熵参数，Fuzdim为模糊熵模式维度
%   -Fuzr：模糊熵参数，Fuzr为模糊熵的阈值
%   -Fuzn：模糊熵参数，Fuzn为模糊熵权重
%   -Pedim：排列熵参数，Pedim为排列熵模式维度
%   -Pet：排列熵参数，Pet为排列熵的时间延迟
%   -fs：采样频率，采样频率即每秒钟采集的数据点数，按照实际情况设置，该参数目前在包络熵特征采集中用到
% featureNamesCell：拟进行特征提取的特征名称，该变量为cell类型，其中包含的特征名称为字符串，特征名称需要在下边列表中：
% 目前支持的特征（2022.7.10，共8种）：
%      psdE：功率谱熵
%      svdpE：奇异谱熵
%      eE：能量熵
%      ApEn：近似熵
%      SampleEn：样本熵
%      FuzzyEn：模糊熵
%      PerEn：排列熵
%      enveEn：包络熵
% 
% 输出：
% fea：数据data的特征值数组，其特征值顺序与featureNamesCell一一对应
% 

%  Copyright (c) 2022 Mr.看海 All rights reserved.
%  代码原理介绍：https://www.zhihu.com/column/xinhao
%  本代码为淘宝买家专用，不开源，请勿公开分享~
rng('default')
data = data';
[len,num] = size(data); %获取数据尺寸，数据长度/数据组数
if len == 1  %一维数据行列方向如果有误
    data = data';
    [len,num] = size(data); %获取数据尺寸，数据长度/数据组数
end
allFeaNames = {'psdE','svdpE','eE','ApEn', 'SpEn','FuzzyEn','PeEn','enveEn'}; %所有特征名称

%% 1.信息熵
% 1.1 功率谱熵
psdE=zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'psdE'}))
    psdE = kPowerSpectrumEntropy(data);
end
% 1.2 奇异谱熵
svdpE =zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'svdpE'}))
    if ~exist('options.svdpEn') %如果未输入窗口值，按照数据长度的1/2取值
        options.svdpEn = round(0.5*(len));
    end
    svdpE = kSingularSpectrumEntropy(data,options.svdpEn);
end
% 1.3 能量熵
eE = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'eE'}))
    eE = kEnergyEntropy(data);
end
allTimeFea = [psdE;svdpE;eE];
%% 2.近似熵指标
ApEn = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'ApEn'}))
     if num == 1
         ApEn = kApproximateEntropy(data, options.Apdim, options.Apr);  %一维数据
    else
        for i = 1:num
            ApEn(i) = kApproximateEntropy(data(:,i), options.Apdim, options.Apr);  %二维数据
        end
    end

end
allTimeFea = [allTimeFea;ApEn];

%% 3.样本熵指标
SpEn = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'SpEn'}))
     if num == 1
         SpEn = kSampleEn(data, options.Spdim, options.Spr);  %一维数据
    else
        for i = 1:num
            SpEn(i) = kSampleEn(data(:,i), options.Spdim, options.Spr);    %二维数据
        end
    end

end
allTimeFea = [allTimeFea;SpEn];

%% 4.模糊熵指标
FuzzyEn = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'FuzzyEn'}))
    if num == 1
        FuzzyEn = kFuzzyEntropy(data,options.Fuzdim,options.Fuzr,options.Fuzn);  %一维数据
    else
        for i = 1:num
            FuzzyEn(i) = kFuzzyEntropy(data(:,i),options.Fuzdim,options.Fuzr,options.Fuzn);  %二维数据
        end
    end
end
allTimeFea = [allTimeFea;FuzzyEn];

%% 5.排列熵指标
PeEn = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'PeEn'}))
     if num == 1
         PeEn = kPermutationEntropy(data, options.Pedim, options.Pet);  %一维数据
     else
        for i = 1:num
            PeEn(i) = kPermutationEntropy(data(:,i), options.Pedim, options.Pet);  %二维数据
        end
     end
end
allTimeFea = [allTimeFea;PeEn];
%% 6.包络熵指标
enveEn = zeros(1,num);%初始化
if sum(contains(featureNamesCell,{'enveEn'}))
     if num == 1
         enveEn = kEnvelopeEntropy(data,options.fs);  %一维数据
     else
        for i = 1:num
            enveEn(i) = kEnvelopeEntropy(data(:,i),options.fs);  %二维数据
        end
     end
end
allFea = [allTimeFea;enveEn];
%% end.对输出变量fea赋值
fea = [];
for i = 1:length(featureNamesCell)
    %featureNamesCell{i} == allFeaNames;
    try
    if find(contains(allFeaNames,featureNamesCell{i})) %如果当前特征在特征列表中
        fea = [fea;allFea(find(strcmp(allFeaNames,featureNamesCell{i})),:)];
    end
    catch ME 
    end
end

fea = fea';

end

function ie = kInformationEntopy(sig,SegmentNum)
% 计算信号的信息熵
% 参考《矿用带式输送机托辊远程故障诊断系统》
% 输入：
% sig：输入信号
% SegmentNum：拟分组数，如果不输入，则自动使用斯特格斯（Sturges）经验公式计算。
% 输出：
% ie：信息熵求解结果
[len,num] = size(sig);
if num >= 2
    SigLen = len;  %输入信号长度
    if nargin == 1
        SegmentNum = round(1.87*(SigLen-1)^(2/5));  %斯特格斯（Sturges）的经验公式，求最佳分组数
    end
    CutLen = SigLen/SegmentNum; %每组信号长度
    Ent = [];
    for i = 1:SegmentNum
        Ent = [Ent;sum(sig(round(CutLen*(i-1)+1):round(CutLen*i),:),1)];  %求每组信号的总能量
    end
    pk = Ent./sum(Ent,2); %第k段能量占总能量大小
    ie = -sum(pk.*log(pk)); %信息熵公式
end 
if num == 1
    
    [SigLen,~] = size(sig);  %输入信号长度
    if nargin == 1
        SegmentNum = round(1.87*(SigLen-1)^(2/5));  %斯特格斯（Sturges）的经验公式，求最佳分组数
    end
    CutLen = SigLen/SegmentNum; %每组信号长度
    for i = 1:SegmentNum
        Ent(i) = sum(sig(round(CutLen*(i-1)+1):round(CutLen*i)));  %求每组信号的总能量
    end
    pk = Ent/sum(Ent); %第k段能量占总能量大小
    ie = -sum(pk.*log(pk)); %信息熵公式
end
end
function svdpE = kSingularSpectrumEntropy(data,n)
% 求信号的奇异谱熵
% 参考《矿用带式输送机托辊远程故障诊断系统》
% 输入：
% data：待分析信号
% n：窗口长度，注意：2<=n<=lenght(data)-1
% 输出：
% svdpE：奇异谱熵值
% data = data(:)'; %强制转化为行向量
% if nargin == 1
%     n = round(sqrt(length(data)));  %斯特格斯（Sturges）的经验公式，求最佳分组数
% end
[len,num] = size(data);
m = len - n - 1; %矩阵A的行数
svdpE = []; 
for j = 1:num
    A = []; %初始化
    for i = 1:m
        A = [A;data(i:i+len-m,j)'];
    end

    svdVal = svd(A); %求奇异值

    svdpE(j) = kInformationEntopy(svdVal,length(svdVal));  %奇异谱熵。分组数等于数据长度。
end
end
function psdE = kPowerSpectrumEntropy(data)
% 求信号的功率谱熵
% 功率谱使用周期图法计算得到
% 参考《矿用带式输送机托辊远程故障诊断系统》
% 输入：
% data：待分析信号
% 输出：
% psdE：功率谱熵值
[len,num] = size(data);
psdE = []; 
for j = 1:num
    [pxx] = periodogram(data(:,j)); %周期图法求功率谱
    [len,~] = size(pxx);
    psdE(j) = kInformationEntopy( pxx,len);
end


end
function eE = kEnergyEntropy(data)
% 求信号基于emd分解算法的能量熵
% 参考《面向高铁走行部故障诊断算法的研究与实现》
% 输入：
% data：待分析信号
% 输出：
% eE：能量熵值

% 需要使用MATLAB2018a及更新版本
[len,num] = size(data);
for i = 1:num
    imf = emd(data(:,i));
    imfE = sum(imf.^2,2);
    eE(i) = kInformationEntopy(imfE,length(imfE));  %奇异谱熵。分组数等于数据长度。
end
end


function ApEn = kApproximateEntropy(data, dim, r)
%  计算近似熵(ApEn)，在网上流传的ApproximateEntropy函数基础上修改
%  输入：
%  data  - 待分析数据，需要是一维数据
%  dim - 模式维度
%  r  -  阈值大小，一般选择r=0.1~0.25，（网上流传的ApproximateEntropy需要再乘以data数据的标准差，本程序不需要，只用输入0.1-0.25间的数字即可）
%  输出：
%  ApEn：近似熵值大小
%  (网上流传的ApproximateEntropy还有一个降采样参数tau，此函数中不进行降采样，删去该参数)
%  论文：Pincus S M . Approximate entropy as a measure of system complexity[J]. Proceedings of the National Academy of Sciences ,1991,88(6):2297—2301.
% 理论讲解：https://zhuanlan.zhihu.com/p/494761890/
% 代码改写与注释由 Mr.看海@知乎 完成

data = data(:);  %强制转化数据为列方向
N = length(data); %获取数据长度
phi = zeros(1,2); %初始化参数
r = r*std(data);%将r的值赋值为比率乘以数据标准差

for j = 1:2
    m = dim+j-1;  %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（5）步。两次循环分别赋值维度为m和m+1
    C = zeros(1,N-m+1);    %C值初始化，网上流传的ApproximateEntropy中此值命名为phi，与理论不对应，故修改为C
    dataMat = zeros(m,N-m+1);    %初始化
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i); %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（1）步，遍历所有指定维度下的折线数据
    end    
    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));  %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（2）步，遍历所有X(i)和X(j)的组合，并求出端点坐标差的绝对值
        boolMat = any( (tempMat > r),1);  %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（2）步，如果X(j)的对应端点都在容限范围内，则认为2维特征向量X(i)和X(j)的模式在r下近似，这里先找到有超出容限的组合，在下一步取非，则可以得到容限内所有组合
        C(i) = sum(~boolMat)/(N-m+1);   %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（3）步，此步取非（即~），可以得到所有阈值容限内的数量
    end
    % summing over the counts
    phi(j) = sum(log(C))/(N-m+1);  %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（4）步，网上流传的ApproximateEntropy中此值命名为result，为了与理论保持一致，变量名改为phi
end
ApEn = phi(1)-phi(2);    %此处对应https://zhuanlan.zhihu.com/p/494761890/中，第（6）步
end

function kSampleEnValue = kSampleEn(data, dim, r)
% 计算data序列的样本熵
% 输入：
% data：一维数据
% dim：重构维数，一般选择1或2，优先选择2，一般不取m>2
% r：阈值大小，一般选择r=0.1~0.25，（网上流传的SampEn需要再乘以data数据的标准差，本程序不需要，只用输入0.1-0.25间的数字即可）
% 输出：
% kSampleEnValue：样本熵值大小
% 论文：Richman J S, Moorman J R. Physiological time-series analysis using approximate entropy and sample entropy[J]. American Journal of Physiology-Heart and Circulatory Physiology, 2000.
% 理论讲解：https://zhuanlan.zhihu.com/p/519042619/
% 代码改写与注释由 Mr.看海@知乎 完成

data = data(:);  %强制转化数据为列方向
N = length(data); %获取数据长度
phi = zeros(1,2); %初始化参数
r = r*std(data); %将r的值赋值为比率乘以数据标准差

for j = 1:2
    m = dim+j-1;  %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（5）步。两次循环分别赋值维度为m和m+1
    B = zeros(1,N-m+1);    %B值初始化，按照论文维度为m和m+1时该变量名应该分别为B和A，但是代码简洁起见，统一使用B表示
    dataMat = zeros(m,N-m+1);    %初始化
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i); %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（1）步，遍历所有指定维度下的折线数据
    end    
    for i = 1:N-m
        tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));  %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（2）步，遍历所有X(i)和X(j)的组合，并求出端点坐标差的绝对值，此时包括了i=j的情况，要在往后数2行的位置减掉这种情况的计数
        boolMat = any( (tempMat > r),1);  %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（2）步，如果X(j)的对应端点都在容限范围内，则认为2维特征向量X(i)和X(j)的模式在r下近似，这里先找到有超出容限的组合，在下一步取非，则可以得到容限内所有组合
        B(i) = (sum(~boolMat)-1)/(N-m-1);   %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（3）步，此步取非（即~），可以得到所有阈值容限内的数量，注意分子上减掉的1是要减掉i=j的情况
    end
    phi(j) = sum(B)/(N-m);  %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（4）步，按照理论此处变量名应为Bm和Am，为了简洁统一命名为phi
end
kSampleEnValue = log(phi(1))-log(phi(2));    %此处对应https://zhuanlan.zhihu.com/p/519042619/中，第（6）步
end

function FuzEn = kFuzzyEntropy(data,dim,r,n)
% 计算模糊熵，在网传FuzzyEntropy.m基础上修改并注释
% 输入：
% data：一维数据
% dim：重构维数，一般选择1或2，优先选择2，一般不取m>2
% r：阈值大小，一般选择r=0.1~0.25，（网上流传的SampEn需要再乘以data数据的标准差，本程序不需要，只用输入0.1-0.25间的数字即可）
% n：权重，取较小的整数，如2或者3。推荐设置为2
% 输出：
% FuzEn：模糊熵值大小

% 论文：Chen W , Wang Z , Xie H , et al. Characterization of Surface EMG Signal Based on Fuzzy Entropy[J]. IEEE Transactions on Neural Systems and Rehabilitation Engineering, 2007, 15(2):266-272.
% 理论讲解：https://zhuanlan.zhihu.com/p/519633544/
% 代码改写与注释由 Mr.看海@知乎 完成

data = data(:)';  %强制转化数据为行向量
N = length(data); %数据长度
phi = zeros(1,2); %初始化参数
r = r*std(data); %将r的值赋值为比率乘以数据标准差

for m = dim:dim+1
    count = zeros(N-m+1,1);    %初始化
    dataMat = zeros(N-m+1,m);  %初始化
    
    % 设置数据矩阵，构造成m维的矢量
    for i = 1:N-m+1
        dataMat(i,:) = data(1,i:i+m-1)-mean(data(1,i:i+m-1));  %%此处对应https://zhuanlan.zhihu.com/p/519633544中，第（1）步，序列定义
    end
    % 利用距离计算相似模式数
    for j = 1:N-m+1
        % 计算切比雪夫距离，不包括自匹配情况
        tempmat=repmat(dataMat(j,:),N-m+1,1);
        dist = max(abs(dataMat - tempmat),[],2); %上行与本行对应https://zhuanlan.zhihu.com/p/519633544中，第（2）步，遍历所有X(i)和X(j)的组合，并求出端点坐标差的绝对值的最大值
        D=exp(-(dist.^n)/r);                     %此处对应https://zhuanlan.zhihu.com/p/519633544中，第（3）步，模糊隶属度
        count(j) = (sum(D)-1)/(N-m-1);           %此处对应https://zhuanlan.zhihu.com/p/519633544中，第（4）步，第一个求和符号内的部分，注意分子中减1是为了减掉i=j情况多出的数值
    end
    phi(m-dim+1) = sum(count)/(N-m);             %此处对应https://zhuanlan.zhihu.com/p/519633544中，第（4）步
end
    % 计算得到的模糊熵值
    FuzEn = log(phi(1)/phi(2));                  %此处对应https://zhuanlan.zhihu.com/p/519633544中，第（5）步
end


function pe = kPermutationEntropy(data,m,t)
% 计算排列熵，在网传PermutationEntropy.m基础上修改并注释
% 输入：
% data：一维数据
% dim：重构维数，推荐值为3~7，如果能够接受程序运行时间久，可以将m设置得尽量大一些。
% t：时间延迟，取整数，设置为1的情况较多。
% 输出：
% pe：排列值大小

% 论文：Bandt C，Pompe B. Permutation entropy:a natural complexity measure for time series[J]. Physical Review Letters,2002,88(17):174102.
% 理论讲解：https://zhuanlan.zhihu.com/p/519809254/
% 代码改写与注释由 Mr.看海@知乎 完成

data = data(:);  %强制转化数据为列方向
N = length(data);  %数据长度
permlist = perms(1:m);  %列出m维度下所有可能的排列
[h,~]=size(permlist);   %所有排列的数量，h=m!
c(1:length(permlist))=0;  %参数初始化，c用于储存各种排序的数量

 for j=1:N-t*(m-1)  %开始循环，注意冒号后为循环排列总数
     [~,iv(j,:)]=sort(data(j:t:j+t*(m-1))); %%此处对应https://zhuanlan.zhihu.com/p/519809254/中，第（1）（2）步，排序使用的是知乎文章中第二种理解方式（括号中的那种）
     for jj=1:h
         if (abs(permlist(jj,:)-iv(j,:)))==0  %遍历各种排列组合permlist，将相同排序记录在变量c的对应位置当中，每有一个对应位置计数加一
             c(jj) = c(jj) + 1 ; %遍历各种排列组合permlist，将相同排序记录在变量c的对应位置当中，每有一个对应位置计数加一
         end
     end
 end
c=c(c~=0);  %剔除c中的0量
p = c/sum(c); %%此处对应https://zhuanlan.zhihu.com/p/519809254/中，第（3）步
pe = -sum(p .* log2(p));%此处对应https://zhuanlan.zhihu.com/p/519809254/中，第（4）步
% 归一化
pe=pe/log2(factorial(m));%此处对应https://zhuanlan.zhihu.com/p/519809254/中，第（5）步
end

function enveEn = kEnvelopeEntropy(data,fs)
% 计算包络熵
% 输入：
% data：一维数据
% fs：采样频率，采样频率即每秒钟采集的数据点数，按照实际情况设置
% 输出：
% enveEn：包络熵值

% 论文《自适应最大相关峭度解卷积方法及其在轴承早期故障诊断中的应用》

% 带通滤波
data = data(:);
ba = [fs*1/4,fs*3/8];  %带通滤波带宽，MATLAB推荐值。https://ww2.mathworks.cn/help/signal/ref/envspectrum.html#namevaluepairarguments
b = fir1(50,[ba(1) ba(2)]/(fs/2));  %带通滤波参数，默认50阶
data = bsxfun(@minus,data,mean(data,1));  %去除直流分量
xBandPass = conv2(data,b(:),'same');  % 进行带通滤波

% 计算包络谱
xAn = hilbert(xBandPass);  %希尔伯特变换
xEnv = abs(xAn);           %求模，即包络值

% 计算包络熵
pj = xEnv./sum(xEnv); %包络熵归一化
enveEn = -sum(pj.*log(pj));
end