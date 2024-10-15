function [imf,CenFs]  = kVMD(data,FsOrT, alpha, K, tol)
% 整合版VMD函数，整合了第三方VMD和MATLAB自带工具箱的VMD分解方法
% 默认设置下按照MATLAB自带库进行
% 并且统一输出格式，配合后续进行希尔伯特-黄变换等无缝衔接
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

% 注意：使用两种库所得到的imf结果会存在差异，这是由两个库自身算法差异导致的，属于正常现象。

% 注意：在使用该代码之前，请务必安装工具箱：http://www.khscience.cn/docs/index.php/2020/04/09/1/
if length(FsOrT) == 1
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
%% 版本判断
vtemp = version;
vv = {vtemp(strfind(vtemp,'R')+1:strfind(vtemp,'R')+4),vtemp(strfind(vtemp,'R')+5)}; %获取版本号
if str2num(vv{1}) < 2020
    type = 2;  %MATLAB版本低，使用第三方库
else
    type = 1;
end


if type == 1
    [imf,~,info] = vmd(data,'AbsoluteTolerance',tol,'NumIMF',K,'PenaltyFactor',alpha);
    % 注意：在MATLAB自带vmd函数的输出参数中，有一项名称为res
    % 然而该res与emd等方法的残差有所不同。
    % emd等方法的残差通常代表信号趋势项，而该处的res则更接近于信号重构误差
    % 为了和其他“类EMD”方法保持一致，此处imf分量中加入res
    % 当前imf的第一个分量对应的是“类EMD”方法中的残差项（也就是趋势项）
    % 
    imf = imf';
    CenFs = info.CentralFrequencies*Fs;
elseif type == 2
    tau=0;      % tau     - time-step of the dual ascent ( pick 0 for noise-slack )
    DC=1;       % DC      - true if the first mode is put and kept at DC (0-freq)
    init=1;     % init    - 0 = all omegas start at 0
                %           1 = all omegas start uniformly distributed
                %           2 = all omegas initialized randomly
    [imf,~,CenFsTemp] = oVMD(data, alpha, tau, K, DC, init, tol);
    CenFs = [CenFsTemp(end,2:end),CenFsTemp(1)]*Fs;
    imfs = imf(2:end,:);
    imf = [imfs;imf(1,:)];
                      %为了和其他“类EMD”方法统一，调整imf排列次序。
                      %如果要得到比较理想的imf排列结果，推荐使用MATLAB2020及以上版本软件
end
%
end