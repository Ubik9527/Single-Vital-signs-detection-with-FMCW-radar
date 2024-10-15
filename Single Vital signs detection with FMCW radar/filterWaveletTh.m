function s = filterWaveletTh(data,wname,SORH,lev,tptr,options)
% 使用小波阈值法实现滤波，本文件是为了形式统一方便调用的一层封装，具体实现在kwden.m、kwdencmp.m、kwthresh.m、kthselect.m文件中
% 这4个文件均为经本人调整过后的函数文件，原函数分别为wden,wdencmp,wthresh,thselect，所以关于这几个函数更多介绍可以查看对应的官方帮助文档
% 输入：
% data：待滤波数据
% wname：小波名称，可选范围参考这里：https://ww2.mathworks.cn/help/wavelet/ref/wfilters.html?searchHighlight=wname&s_tid=srchtitle_wname_2#d123e130597
% SORH：阈值函数类型
%       sorh ='a1'时，采用改进方法1，改进后的软阈值：
%                     参考论文：孙万麟, 王超. 基于改进的软阈值小波包网络的电力信号消噪[J]. 海军工程大学学报, 2019(4).
%       sorh = 'a2'时，采用改进方法2：
%                     参考论文：刘冲, 马立修, 潘金凤,等. 联合VMD与改进小波阈值的局放信号去噪[J]. 现代电子技术, 2021, 44(21):6.
%       sorh = 'a3'时，采用改进方法4：
%                     参考论文：基于改进小波阈值-CEEMDAN算法的ECG信号去噪研究
%       sorh = 'a4'时，采用改进方法3：
%                     参考论文：基于改进小波阈值函数的语音增强算法研究
%                     采用该方法需要输入两个调节因子，分别是alpha和gamma，取alpha>0，0<gamma<1
%       sorh = 'a5'时，采用改进方法3：
%                     参考论文：《基于VMD与改进小波阈值的地震信号去噪方法研究》
%                     采用该方法需要输入两个调节因子，分别是alpha和beta% lev：  小波分解水平
% tptr： 阈值选择规则，有可选类型：
%       '' — Adaptive threshold selection using the principle of Stein's Unbiased Risk Estimate (SURE).
%       'sqtwolog' — Fixed-form threshold is sqrt(2*log(length(X))).
%       'heursure' — Heuristic variant of 'rigrsure' and 'sqtwolog'.
%       'minimaxi' — Minimax thresholding.
%       'visushrink' - 通用阈值去噪方法
% options：部分阈值函数需要补充的参数设置，即某些调节因子，需要用结构体方式调用幅值，如options.a4_alpha=2。如果不需要设置，可以赋值为options=[]。具体如下：
%          -a4_alpha：改进方法4的alpha值设置
%          -a4_gamma：改进方法4的gamma值设置
%          -a5_alpha：改进方法5的alpha值设置
%          -a5_beta：改进方法5的beta值设置
% 输出：
% s：经过滤波后得到的数据

%  Copyright (c) 2022 Mr.看海 All rights reserved.
%  本代码为淘宝买家专用，不开源，请勿公开分享~

s = kwden(data,tptr,SORH,'one',lev,wname,options);   %调用了经 Mr.看海 修改后的小波阈值滤波函数
end