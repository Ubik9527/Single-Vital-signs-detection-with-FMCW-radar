function y = kwthresh(x,sorh,t,options)
% �� Mr.���� �޸Ĺ������ֵ������ԭ����Ϊwthresh��
% ����˸Ľ���ֵ��������
% ���룺
% x����ֵ����������ֵ�������޸�
% sorh����ֵ��������
%       sorh ='a1'ʱ�����øĽ�����1���Ľ��������ֵ��
%                     �ο����ģ�������, ����. ���ڸĽ�������ֵС��������ĵ����ź�����[J]. �������̴�ѧѧ��, 2019(4).
%       sorh = 'a2'ʱ�����øĽ�����2��
%                     �ο����ģ�����, ������, �˽��,��. ����VMD��Ľ�С����ֵ�ľַ��ź�ȥ��[J]. �ִ����Ӽ���, 2021, 44(21):6.
%       sorh = 'a3'ʱ�����øĽ�����3��
%                     �ο����ģ����ڸĽ�С����ֵ-CEEMDAN�㷨��ECG�ź�ȥ���о�
%       sorh = 'a4'ʱ�����øĽ�����4��
%                     �ο����ģ����ڸĽ�С����ֵ������������ǿ�㷨�о�
%                     ���ø÷�����Ҫ���������������ӣ��ֱ���alpha��gamma��ȡalpha>0��0<gamma<1
%       sorh = 'a5'ʱ�����øĽ�����5��
%                     �ο����ģ�������VMD��Ľ�С����ֵ�ĵ����ź�ȥ�뷽���о���
%                     ���ø÷�����Ҫ���������������ӣ��ֱ���alpha��beta
% t����ֵ
% �����
% y����ֵ�������ֵ
% ������ֵ����ͼ�߿������������룺
% figure('color','w');plot(-5:0.1:5,kwthresh(-5:0.1:5,'a1',1));
% ������ֵ����ͼ������Ӳ��ֵͼ�߶Ա�ͼ�����������룺
% figure('color','w');plot(-5:0.1:5,kwthresh(-5:0.1:5,'a1',1));hold on
% plot(-5:0.1:5,kwthresh(-5:0.1:5,'s',1))
% plot(-5:0.1:5,kwthresh(-5:0.1:5,'h',1))
% legend('�Ľ���ֵ����','����ֵ','Ӳ��ֵ')

%  ������Ϊ�Ա����ר�ã�����Դ�����𹫿�����~

%WTHRESH Perform soft or hard thresholding. 
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's')
%   or hard (if SORH = 'h') T-thresholding  of the input 
%   vector or matrix X. T is the threshold value.
%
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft 
%   thresholding is shrinkage.
%
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard
%   thresholding is cruder.
%
%   See also WDEN, WDENCMP, WPDENCMP.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 24-Jul-2007.
%   Copyright 1995-2020 The MathWorks, Inc.

if isStringScalar(sorh)
    sorh = convertStringsToChars(sorh);
end

switch sorh
  case 's'
    tmp = (abs(x)-t);
    tmp = (tmp+abs(tmp))/2;
    y   = sign(x).*tmp;
 
  case 'h'
    y   = x.*(abs(x)>t);
  case 'a1'
        y = (abs(x)>t).*sign(x).*sqrt(abs(x).^2-t^2);
  case 'a2'
    for i = 1:length(x)
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))-2^(t-abs(x(i))));
        else
            y(i) = 0;
        end
    end
  case 'a3'
    for i = 1:length(x)
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))-2*t/(exp((abs(x(i))-t)/t)+1));
        else
            y(i) = 0;
        end
    end
  case 'a4'  %�ο� ���ڸĽ�С����ֵ������������ǿ�㷨�о�
        alpha = options.a4_alpha; %��������
        gamma = options.a4_gamma; %��������
        for i = 1:length(x)
            if abs(x(i)) > t
                y(i) = sign(x(i))*(abs(x(i))-(t/(1+alpha).*gamma.^(sqrt(x(i).^2-t.^2))));
            else
                y(i) = sign(x(i))*alpha/(1+alpha)*exp(10*(abs(x(i))-t))*abs(x(i));
            end
        end

  case 'a5' %�ο� ������VMD��Ľ�С����ֵ�ĵ����ź�ȥ�뷽���о���
    for i = 1:length(x)
        alpha = options.a5_alpha; %��������
        beta = options.a5_beta; %��������
        if abs(x(i)) > t
            y(i) = sign(x(i))*(abs(x(i))- t/power((abs(x(i)).^beta-(abs(t)).^beta+1),1/beta)*1/exp(power(abs(x(i))-abs(t),1/alpha)));
        else
            y(i) = 0;
        end
    end
    
  otherwise
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end
