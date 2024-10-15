function reSig = filEMDsWaveletTh(sig)
% ICEEMDAN-������-�Ľ���С����ֵ�ķ�װ�������ú����е�ICEEMDAN�����������ء�С����ֵ����ز������Խ����޸��滻
% �������:
%   sig - ���˲��ź�
%
% �������:
%   reSig - �˲�����ź�

%% 1.ICEEMDAN�ֽ�
% Nstd = 0.2; %NstdΪ����������׼����Y��׼��֮��
% NE = 100;   %NEΪ���źŵ�ƽ������
% MaxIter = 1000;% MaxIter ����������
% imf = pICEEMDAN(sig,1,Nstd,NE,MaxIter);
FsOrT = 50;
alpha = 2000;
K = 6;
tol = 1e-6;
imf = pVMDandFFT(sig,FsOrT, alpha, K, tol);
%% 2.��������ȡ�����Ҫ�������������������Բο�testGenFeaEn.m�ļ�����ʾ
featureNamesCell = {'SpEn'}; %Ҫ����������ȡ���������ƣ����Ը���ʵ����Ҫ����ɾ�������µ�����ע��ƴд��ȷ

% ���������ز�����SpdimΪ�����ص�ģʽά�ȣ�SprΪ�����ص���ֵ���������ȡ��������������ɾ����������
option.Spdim  = 2;
option.Spr   = 0.15;

fea = genFeatureEn(imf,featureNamesCell,option);  %����genFeature���������������ȡ�����������ֵ�ᱣ����fea�����
                                             %fea�����ĳ��Ⱥ�featureNamesCell��ָ����������һ�£���˳��һһ��Ӧ
                                             %����������ɺ���MATLAB�Ĺ�������˫��fea���������Բ鿴��õľ�����ֵ
% disp(['������ֵΪ:',num2str(fea')])
% % ��������ֵͼ�����ڷ�����IMF�����ĸ�����
% figure('Color','w');
% bar(fea);  % ��������ͼ��ʾÿ��IMF�����Ķ�߶�������
% hold on;
% plot(fea, 'r--o');  % ������ͼ�ϸ��ǻ�������ͼ���Ա�ֱ�۱Ƚ�
% title(featureNamesCell{1});
% xlabel('IMF');
% ylabel('��ֵ');

% %% 3.������ֵɸѡIMF����
% th = 0.1; %��ֵ
% selectedIMF = imf(fea>th,:); % ɸѡ�����ش���0.1��IMF����
% selectedIdx = find(fea>th); % �ҳ�����������IMF�������
% selectedIMF = imf(1:5,:);
% selectedIdx = [1,2,3,4,5].';
% % fprintf('����������IMF�������Ϊ:');
% % disp(selectedIdx');

%% 3. ����IMF��ԭ�źŵ������ɸѡIMF����
thre = 0.5; %��ֵ
for i = 1:K
    R = corrcoef(sig,imf(i,:));
    r(i) = R(1,2);
end
fprintf('IMF���������ϵ��:');
disp(r);
selectedIMF = imf(r<thre,:); % ɸѡ�����ش���0.1��IMF����
selectedIdx = find(r<thre); % �ҳ�����������IMF�������
fprintf('��Ҫ�˲���IMF����Ϊ:');
disp(selectedIdx);
%% 3.����Ƶ��ɸѡIMF����
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
% fprintf('����������IMF�������Ϊ:');
% disp(selectedIdx);
%% 4.��ѡ����IMF��������С����ֵ�˲�
wname = 'sym10'; % С������
SORH = 'a3'; % ��ֵ����,���øĽ�������ֵa4
lev = 4; % �ֽ����
tptr = 'heursure'; % ��ֵѡ�����
options.a4_alpha = 2; % a4������alpha����
options.a4_gamma = 0.8; % a4������gamma����

denoisedIMFs = zeros(size(selectedIMF)); % ��ʼ��

N_selected = size(selectedIMF,1);
figure('Color','w');
N_selected = size(selectedIMF,1);
for i=1:N_selected
    denoisedIMFs(i,:) = filterWaveletTh(selectedIMF(i,:),wname,SORH,lev,tptr,options);
    % �����˲�ǰ��Ա�ͼ
    subplot(N_selected,2,2*i-1);
    plot(selectedIMF(i,:),'k');
    title(['IMF ',num2str(selectedIdx(i)),' (�˲�ǰ)']);
    ylim_curr = ylim; % ��¼��ǰ��ͼ��y�᷶Χ
    
    subplot(N_selected,2,2*i);
    plot(denoisedIMFs(i,:),'k');
    title(['IMF ',num2str(selectedIdx(i)),' (�˲���)']);
    ylim(ylim_curr); % ��y�᷶Χ����Ϊ�������ͼ��ͬ
end
%% 5.�ع��˲�����ź�
reSig = sum(denoisedIMFs,1);

% ��δɸѡ��IMF����ֱ����ӵ������
unselectedIMF = imf(~ismember(1:size(imf,1),selectedIdx),:);
reSig = reSig + sum(unselectedIMF,1);

end