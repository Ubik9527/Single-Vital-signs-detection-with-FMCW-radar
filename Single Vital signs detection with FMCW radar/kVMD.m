function [imf,CenFs]  = kVMD(data,FsOrT, alpha, K, tol)
% ���ϰ�VMD�����������˵�����VMD��MATLAB�Դ��������VMD�ֽⷽ��
% Ĭ�������°���MATLAB�Դ������
% ����ͳһ�����ʽ����Ϻ�������ϣ������-�Ʊ任���޷��ν�
% ���룺
% data�����ֽ������(һά)
% FsOrT������Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά���������δ֪����Ƶ�ʣ�������Ϊ1
% alpha   - �ͷ�����
% K       - ָ���ֽ�ģ̬��
% tol     - �����ݲ���Ż���ֹͣ׼��֮һ������ȡ 1e-6~5e-6

% �����
% imf���ں�ģ̬������ͳһΪn*m��ʽ������nΪģ̬����mΪ���ݵ��������� imf(1,:)��IMF1��imf(end,:)��Ϊ�в�
% ע�⣬Ϊ������������EMD�������ֽ������imf��������һ�£��������ڽ�imf��������˷�ת
% ����֤imf�����ǴӸ�Ƶ���Ƶ����
% CenFs����CentralFrequencies����imf����������Ƶ��

% ע�⣺ʹ�����ֿ����õ���imf�������ڲ��죬�����������������㷨���쵼�µģ�������������

% ע�⣺��ʹ�øô���֮ǰ������ذ�װ�����䣺http://www.khscience.cn/docs/index.php/2020/04/09/1/
if length(FsOrT) == 1
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
%% �汾�ж�
vtemp = version;
vv = {vtemp(strfind(vtemp,'R')+1:strfind(vtemp,'R')+4),vtemp(strfind(vtemp,'R')+5)}; %��ȡ�汾��
if str2num(vv{1}) < 2020
    type = 2;  %MATLAB�汾�ͣ�ʹ�õ�������
else
    type = 1;
end


if type == 1
    [imf,~,info] = vmd(data,'AbsoluteTolerance',tol,'NumIMF',K,'PenaltyFactor',alpha);
    % ע�⣺��MATLAB�Դ�vmd��������������У���һ������Ϊres
    % Ȼ����res��emd�ȷ����Ĳв�������ͬ��
    % emd�ȷ����Ĳв�ͨ�������ź���������ô���res����ӽ����ź��ع����
    % Ϊ�˺���������EMD����������һ�£��˴�imf�����м���res
    % ��ǰimf�ĵ�һ��������Ӧ���ǡ���EMD�������еĲв��Ҳ���������
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
                      %Ϊ�˺���������EMD������ͳһ������imf���д���
                      %���Ҫ�õ��Ƚ������imf���н�����Ƽ�ʹ��MATLAB2020�����ϰ汾���
end
%
end