function [imf,CenFs] = pVMD(y,FsOrT, alpha, K, tol)
% ���ź�VMD�ֽ�ͼ
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
%% 1.���ƶ�άͼ
figure('Name','VMD�ֽ����IMF��������ͼ','Color','white');
subplot(size(imf,1)+1,1,1);
plot(t,y,'k');grid on;
ylabel('ԭʼ����');
title('VMD�ֽ�');
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

%% 2.������άͼ
figure('Name','ģ̬�ֽ����IMF��������ͼ����ά��','Color','white');
imfall = [y(:)'; imf];
[m, n] = size(imfall);
[X, Y] = meshgrid(1 : m, 1 : n );    % ����xyƽ����������X��Y����n��m�еľ���X��ÿһ�ж���(1:m)��Y��ÿһ�ж���(n:-1:1)
Z = imfall(1 : m, :);                    % ѡ��Ҫ��ʾ��ʱ���ź����ݣ�Z��m��n�еľ���

plot3(X, Y, Z); grid on;

set(gca,'Ydir','reverse'); % ��תy��
xticks(1:m)
xticklabels(['ԭʼ�ź�', strcat('IMF', string(1:m-2)), 'res']); % ����x��Ŀ̶ȱ�ǩ

ylabel('��������');
% ylabel('ʱ��(s)');
% zlabel('��ֵ(W)');
zlabel('��ֵ');
set(gcf,'Color', 'white');

title('VMD�ֽ�-��άչ��ͼ');

view(-35, 30);
end
