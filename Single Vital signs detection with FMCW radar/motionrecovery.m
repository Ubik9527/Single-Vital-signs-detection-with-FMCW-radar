function y = motionrecovery(x,thre)

%% 二阶差分计算  
acc = diff(x,2);   %取数组的两阶差分作为加速度
acc = padarray(acc,[0 1],'replicate','both'); %填充因为差分而少掉的数组

%% 人体运动趋势指标
% 取n秒作为一个窗口，将窗口内的所有加速度值取均值，用该值判断窗口内是否发生了移动
fa_divtime = 50*1;    
fa_winnun = 3000/fa_divtime;
fa = zeros(1,fa_winnun);
for k=1:fa_winnun
    for i=1:fa_divtime
        fa(k) = fa(k) + abs(acc((k-1)*fa_divtime+i))/fa_divtime;
        fa_var(k) = var(acc((k-1)*fa_divtime+1:k*fa_divtime));
    end 
end
fa_average = mean(fa);
fa_var_average = mean(fa_var);
figure;
X = 0.02:0.02:60;
plot(X,acc);hold on;
plot(fa,'r');hold on
plot([0,60],[0.2,0.2],'r');
% hold on;plot(acc_divtime:acc_divtime:3000,v);
xlabel('时间(s)');
ylabel('幅值')
title('人体运动检测');
legend('二阶差分值','运动趋势函数值');
% X = 1:acc_divtime:3000;

% figure;
% plot(fa);hold on
% plot([0,60],[0.2,0.2],'r');hold on;
% plot([0.60],[0.01,0.01],'r');
% title('人体运动趋势函数');
% %% 二阶差分失真指标
% u = acc;
% thre=0.5;alpha=-100;
% for k=1:length(angle_fft_last)
%     fu(k)=exp(alpha*(u(k)-thre));
%     (1+exp(alpha*(u(k)-thre)));
% end
% figure;
% plot(fu);
% title('二阶差分失真补偿指标');


%% 判断哪里发生了运动
    %阈值取二阶差分数据的平均数，大于该值的被判断为发生运动部分
    move_index = [];
    indicesToDelete = [];
    for i=10:fa_winnun
        if(fa(i) > thre || fa(i) < 0.01  )
            move_index(end+1) = i; 
        end
    end

    if(isempty(move_index))
        disp('人体无运动发生');
    else
        move_newindex = move_index;
    %     防止误判，将1s以内的小幅运动去除
        for i=2:length(move_index)-1
            if(move_index(i+1)-move_index(i))>1 && (move_index(i)-move_index(i-1))>1
%                 move_newindex(i) = []; 
                indicesToDelete(end+1) = i; 
            end
        end

    
        N = length(move_newindex);
        if(N>2)
            if(move_newindex(N)-move_newindex(N-1)>1)
                indicesToDelete(end+1) = N; 
            end
            if(move_newindex(2)-move_newindex(1)>1)
                indicesToDelete(end+1) = 1; 
            end
        else 
            indicesToDelete(end+1) = N;
            indicesToDelete(end+1) = 1;
        end
        % 创建一个逻辑索引向量，标记要删除的元素的位置
        logicalIndices = false(size(move_newindex));
        logicalIndices(indicesToDelete) = true;
        
        % 使用逻辑索引删除元素
        move_newindex(logicalIndices) = [];
        disp(['检测到运动部分索引',num2str(move_index)]);
        disp(['去除单独索引后',num2str(move_newindex)]);
        
        if(isempty(move_newindex))
            disp('人体无运动发生');
        else
            %将运动分段，并标记
            % tstart1 = (move_newindex(1)-1)*acc_divtime;
            tse_index = [];
            tse_index(end+1) = move_newindex(1);
            for(i=2:length(move_newindex))
                if(move_newindex(i)-move_newindex(i-1)>1)
                    tse_index(end+1) = move_newindex(i-1);
                    tse_index(end+1) = move_newindex(i);
                end
            end
            tse_index(end+1) = move_newindex(end);
            
            for i=1:length(tse_index)/2
                tstart(i) = (tse_index((i-1)*2+1)-1)*fa_divtime;
                tend(i) = (tse_index(i*2)+1)*fa_divtime;
                if(tend(i)>3000)
                    tend(i)=tend(i)-fa_divtime;
                end
                disp(['运动开始',num2str(tstart(i)),'运动结束:',num2str(tend(i))]);
            end
    
            %% 粗略判断信号周期
            angle_fft_rough = diff(x,1);  %相位差分
            angle_fft_rough = padarray(angle_fft_rough,[0 1],'replicate','pre'); %补足数据
            for i=2:length(angle_fft_rough)  %相位突变噪声去除
                if(abs(angle_fft_rough(i))>5*mean(abs(angle_fft_rough)))
                    angle_fft_rough(i)=angle_fft_rough(i-1);
                end
            end
            angle_fft_rough=smoothdata(angle_fft_rough,'movmean',5);  %数据平滑
            N_rough=length(angle_fft_rough);fs = 50;
            FFT_rough = abs(fft(angle_fft_rough));            %--FFT   取模，幅度
            f=(0:N_rough-1)*(fs/N_rough);             %其中每点的频率
            %傅里叶变换结果对称
            figure;
            plot(f(1:N_rough/8),FFT_rough(1:N_rough/8)) %取前一部分放大观察
            xlabel('频率（Hz）');
            ylabel('幅度');
            title('相位信号FFT  ');
            breathmax_rough = 0;
            [~, index_rough] = max(FFT_rough(8:30));
            index_rough = index_rough+7;
            period = round(N_rough/(index_rough-1));

            %% SARIMA补偿
            %启动计时器
            tic
            % 1.加载数据
            data_raw = x.';
            for(i=1:length(tse_index)/2)
               
                step = tend(i)-tstart(i); 
                data = data_raw(1:tstart(i)); 
                S = period; %季节性序列变化周期
                while(length(data)<4*step||length(data)<5*S)
                    data = padarray(data,S,'symmetric','pre'); %当训练数据过短无法拟合模型时，复制一个周期的数据在前面
                end
            
            
                % 2.确定季节性与非季节性差分数，D取默认值1，d从0至3循环，平稳后停止
                dY = data;
                d=0;
                while(adftest(dY)==0) %通过差分次数来判断 
                    dY = diff(dY);
                    d = d+1;
                end
                % 3.确定阶数ARlags,MALags,SARLags,SMALags
                
                AR_Order = 1;
                MA_Order = 1;
                SAR_Order = 1;
                SMA_Order = 1;
                Mdl = creatSARIMA(AR_Order,MA_Order,SAR_Order,SMA_Order,S,d);  %创建SARIMA模型 
                try
                    EstMdl = estimate(Mdl,data);
                catch ME %捕捉错误信息
                    msgtext = ME.message;
                    if (strcmp(ME.identifier,'econ:arima:estimate:InvalidVarianceModel'))
                         msgtext = [msgtext,'  ','无法进行arima模型估计，这可能是由于用于训练的数据长度较小，而要进行拟合的阶数较高导致的，请尝试减小max_ar和max_ma的值']
                    end
                    msgbox(msgtext, '错误')
                    return
                end
                %% 残差检验
                [res,~,logL] = infer(EstMdl,data);   %res即残差
                
                stdr = res/sqrt(EstMdl.Variance);
%                 figure('Name','残差检验','Visible','on')
%                 subplot(2,2,1)
%                 plot(stdr)
%                 title('标准化残差')
%                 subplot(2,2,2)
%                   qqplot(stdr)
%                 title('残差分布')
%                 subplot(2,2,3)
%                 autocorr(stdr)
%                 title('残差自相关函数')
%                 subplot(2,2,4)
%                 parcorr(stdr)
%                 title('残差偏自相关函数')
            
            
                % Durbin-Watson 统计是计量经济学分析中最常用的自相关度量
                diffRes0 = diff(res);  
                SSE0 = res'*res;
                DW0 = (diffRes0'*diffRes0)/SSE0; % Durbin-Watson statistic，该值接近2，则可以认为序列不存在一阶相关性。
                % 5.预测
                [forData,YMSE] = forecast(EstMdl,step,data);   %matlab2018及以下版本写为Predict_Y = forecast(EstMdl,step,'Y0',Y);   matlab2019写为Predict_Y = forecast(EstMdl,step,Y);
                
%                 lower = forData - 1.96*sqrt(YMSE); %95置信区间下限
%                 upper = forData + 1.96*sqrt(YMSE); %95置信区间上限          
%                 figure('Visible','on')
%                 plot(data,'Color',[.7,.7,.7]);
%                 hold on
%                 h1 = plot(length(data):length(data)+step,[data(end);lower],'r:','LineWidth',2);
%                 plot(length(data):length(data)+step,[data(end);upper],'r:','LineWidth',2)
%                 h2 = plot(tstart(i):tend(i),[data(end);forData],'k','LineWidth',2);
%                 plot(data_raw);
%                 legend([h1 h2],'95% 置信区间','预测值','原信号','Location','NorthWest')
%                 title('Forecast')
%                 hold off
                %对预测数据最后一个与原数据相差过大处理
                dif = data_raw(tend(i)) - forData(end); %连续值之间的相位差
                data_raw(tend(i):end) = data_raw(tend(i):end) - dif;
                 %数据替换
                data_raw(tstart(i)+1:tend(i)) = forData(1:end);
            end 
            elapsed_time = toc;
            disp(['SARIMA模型修复运行时间：', num2str(elapsed_time), '秒']);
            figure;
            X = 0.02:0.02:60;
            plot(X,x);
            hold on;plot(X,data_raw,'r');
            title('运动失真修复前后对比');
            legend('原始生命体征信号','恢复后生命体征信号');
            xlabel('时间(s)');ylabel('相位(rad)');
            x = data_raw;
        end
    end
    y = x;
end
