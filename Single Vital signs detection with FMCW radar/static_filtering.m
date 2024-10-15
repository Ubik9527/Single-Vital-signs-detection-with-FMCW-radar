%% 对距离FFT数据进行静态滤波
% 输入参数：
% data_input：距离FFT后的数据
% rangeRes：距离分辨率
% n：选择静态滤波方法
%   0：不滤波
%   1：脉冲对消法
%   2：平均相消法
%   3：相位判断法
% 输出参数：
% data_out：静态滤波后数据
function data_out = static_filtering(data_input,n,rangeRes)
    switch n
        case 1
        %% 1.脉冲对消法
            [rangelen,chriplen] = size(data_input);
            datatmp2 = zeros(rangelen,chriplen-1);
            datatmp2 = complex(datatmp2,datatmp2);
            for j = 1:chriplen-1
                datatmp2(:,j) = data_input(:,j+1) - data_input(:,j);
            end
            % 最后一列用最后一列的数据减去第一列的数据？（或者直接不用）
            DataStaticRemoved_pulseoffset = padarray(datatmp2,[0 1],'replicate','post');
            fft_data_pulseoffset = abs(DataStaticRemoved_pulseoffset).';
            
            % 三维图生成
            fft1d_pulseoffset= zeros(chriplen,rangelen);
                for b=1:chriplen
                    fft1d_pulseoffset(b,:) = (fft_data_pulseoffset(b,:));%
                end
            [X,Y] = meshgrid((0:rangelen-1)*rangeRes, ...
               (0.02:0.02:60)); 
            % 原始数据
            fftdata_raw = abs(data_input);
            fftdata_raw = fftdata_raw.';
            %滤波前后比较
            figure;
            subplot(2,1,1);
            mesh(X,Y,fftdata_raw);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT before static filtering');xlim([0 3]);
            subplot(2,1,2);
            mesh(X,Y,fft1d_pulseoffset);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT after static filtering');xlim([0 3]);
            data_out = DataStaticRemoved_pulseoffset;
        case 2
        %% 平均相消算法
            [rangelen,chriplen] = size(data_input);
            datatmp2 = data_input;
            for i = 1:rangelen
                meandata = mean(datatmp2(i,:));
                datatmp2(i,:) = datatmp2(i,:) - meandata;
            end
            DataStaticRemoved_meanoffset = datatmp2;
            fft_data_meanoffset = abs(DataStaticRemoved_meanoffset).';
            
            % 三维图生成
            fft1d_meanoffset= zeros(chriplen,rangelen);
                for b=1:chriplen
                    fft1d_meanoffset(b,:) = (fft_data_meanoffset(b,:));%
                end
            [X,Y] = meshgrid((0:rangelen-1)*rangeRes, ...
               (0.02:0.02:60));   
            % 原始数据
            fftdata_raw = abs(data_input);
            fftdata_raw = fftdata_raw.';
            %滤波前后比较
            figure;
            subplot(2,1,1);
            mesh(X,Y,fftdata_raw);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT before static filtering');xlim([0 3]);
            subplot(2,1,2);
            mesh(X,Y,fft1d_meanoffset);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT after static filtering');xlim([0 3]);
            data_out = DataStaticRemoved_meanoffset;
        case 3
        %% 静态滤波：找到距离的几个极值点并去掉其中静止的背景
            % 原始数据
            fftdata_raw = abs(data_input);
            fftdata_raw = fftdata_raw.';
            % 求最大距离门
            range_profile_abs=abs(data_input);
            [N,M] = size(data_input);
            range_profile_extreme = zeros(N,1);
            range_profile_extreme_data = zeros(N,M);
            range_profile_extreme_data = complex(range_profile_extreme_data,range_profile_extreme_data);
            range_profile_last=zeros(N,1);
            range_max=1;range_second=1;
            for j = 1:N
               if((j*rangeRes)<2.5 &&(j*rangeRes)>0.5) % 限定了检测距离为0.5-2.5m
                    for i = 1:M             % 进行非相干积累
                        range_profile_last(j) = range_profile_last(j) + power(range_profile_abs(j,i),2);
                    end
                    if (range_profile_last(j)>range_profile_last(range_max))
                        range_second=range_max;
                        range_max=j;
                    elseif (range_profile_last(j)>range_profile_last(range_second))
                        range_second=j;
                    end
                end
            end 
            %取能量最大的几个极值点所在位置组成备选数组
            range_profile_count = 0; %数组下标
            for i=3:N-1
                if range_profile_last(i) > range_profile_last(i-1) && range_profile_last(i )> range_profile_last(i+1) 
                    range_profile_count = range_profile_count+1;
                    range_profile_extreme(range_profile_count) = i;
                end
            end
            range_profile_var = zeros(range_profile_count,1);
            %计算这些点的方差
            for i=1:range_profile_count
                range_profile_extreme_data(i,:) = data_input(range_profile_extreme(i),:);
                range_profile_real = real(range_profile_extreme_data);
                range_profile_imag = imag(range_profile_extreme_data);
                range_profile_phase = atan2(range_profile_imag,range_profile_real);
                range_profile_var = var(range_profile_phase,0,2);
            end
            
            %方差小于一个阈值则视为静止物体,暂定0.5,则在距离FFT上将数据置零，即除去该物体
            for i=1:range_profile_count
                if(range_profile_var(i)<1)
                   data_input(range_profile_extreme(i),:)=0;
                end
            end
            % 去除直流分量
            data_input(1:4,:) = 0;  %置零，去除直流分量
            DataStaticRemoved_variancejudge = data_input;
            fft_data_variancejudge = abs(data_input(:,:));
            fft_data_variancejudge = fft_data_variancejudge.';%非共轭翻转1024*256
            % fft_data_new(:,1:4)=0;         %置零，去除直流分量
            
            % 三维图生成
            fft1d_variancejudge= zeros(M,N);
                for b=1:M
                    fft1d_variancejudge(b,:) = (fft_data_variancejudge(b,:));%
                end
            [X,Y] = meshgrid((0:N-1)*rangeRes, ...
               (0.02:0.02:60));   
            figure;
            subplot(2,1,1);
            mesh(X,Y,fftdata_raw);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT before static filtering');xlim([0 3]);
            subplot(2,1,2);
            mesh(X,Y,fft1d_variancejudge);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT after static filtering');xlim([0 3]);
            data_out = DataStaticRemoved_variancejudge;
        otherwise
            [rangelen,chriplen] = size(data_input);
            fftdata_raw = abs(data_input);
            fftdata_raw = fftdata_raw.';
            [X,Y] = meshgrid((0:rangelen-1)*rangeRes,(0.02:0.02:60)); 
            figure;
            mesh(X,Y,fftdata_raw);view(2);
            xlabel('range(m)');ylabel('time(s)');zlabel('amplitude');
            title('1D-FFT before static filtering');xlim([0 3]);
            data_out = data_input;
    end
    
end