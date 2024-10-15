function [imf,CenFs] = pVMDandFFT(y,FsOrT, alpha, K, tol)
% ���ź�VMD�ֽ����IMF����Ƶ�׶���ͼ
% ���룺
% y�����ֽ������(һά)
% FsOrT�� ����Ƶ�ʻ����ʱ�����������Ϊ����Ƶ�ʣ��ñ������뵥��ֵ�����Ϊʱ���������ñ���Ϊ��y��ͬ���ȵ�һά���������δ֪����Ƶ�ʣ�������Ϊ1
% alpha   - �ͷ�����
% K       - ָ���ֽ�ģ̬��
% tol     - �����ݲ���Ż���ֹͣ׼��֮һ������ȡ 1e-6~5e-6

% �����
% imf���ں�ģ̬������ͳһΪn*m��ʽ������nΪģ̬����mΪ���ݵ��������� imf(1,:)��IMF1��imf(end,:)��Ϊ�в�
% ע�⣬Ϊ������������EMD�������ֽ������imf��������һ�£��������ڽ�imf��������˷�ת
% ����֤imf�����ǴӸ�Ƶ���Ƶ����
% CenFs����CentralFrequencies����imf����������Ƶ��

% ��1����FsOrTΪ����Ƶ�ʣ�
% fs = 100;
% t = 1/fs:1/fs:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pVMDandFFT(y,fs,2000,2,1e-7);
% ��2����FsOrTΪʱ����������Ҫע���ʱFsOrT�ĳ���Ҫ��y��ͬ��
% t = 0:0.01:1;
% y = sin(2*pi*5*t)+2*sin(2*pi*20*t);
% imf = pVMDandFFT(y,t,2000,2,1e-7);

% ע�⣺��ʹ�øô���֮ǰ������ذ�װ�����䣺http://www.khscience.cn/docs/index.php/2020/04/09/1/

%  Copyright (c) 2021 Mr.���� All rights reserved.
%  ������Ϊ�Ա����ר�ã�����Դ�����𹫿�����~
%%
if length(FsOrT) == 1
    t = 1/FsOrT:1/FsOrT:length(y)/FsOrT;
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
[imf,CenFs] = kVMD(y, FsOrT,alpha, K, tol);
figure('Name','VMD�ֽ����IMF����Ƶ�׶���ͼ','Color','white');
subplot(size(imf,1)+1,2,1);
plot(t,y,'k');grid on;
ylabel('ԭʼ����');
title('VMD�ֽ�');
set(gca,'XTick',[]);
subplot(size(imf,1)+1,2,2);
pFFT(y,Fs);grid on;
title('��ӦƵ��');
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
