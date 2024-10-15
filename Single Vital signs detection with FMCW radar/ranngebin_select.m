%% 不同距离门选择函数提取生命体征信号
% 输入参数：
% data_out：距离FFT后信号
% fs：信号采样频率
% t：分段距离门时间间隔单位秒（3,4方法用得着）
% rangeRes：距离分辨率
% n：距离门选择方法：
%     1：固定距离门
%     2：最大距离门
%     3：根据一段时间内的距离门单元众数选择
%     4：根据一段时间内的非相干累积值选择
% 输出参数：
% vital_data：提取到的生命体征信号
% range_index：选中距离门索引

function [vital_data,range_index] = ranngebin_select(data_out,fs,t,rangeRes,n)
    [N,M] = size(data_out);
    real_data = real(data_out);%实部
    imag_data = imag(data_out);%虚部
    
    for i = 1:N
        for j = 1:M  %对每一个range bin取相位 extract phase（弧度rad）
            angle_fft(i,j) = atan2(imag_data(i, j),real_data(i, j));
        end
    end
    angle_fft = angle_fft.';
    switch n
        case 1 
        %% 1.固定距离门
            % 求最大距离门
            range_profile_abs=abs(data_out);
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
            angle_fft_last = angle_fft(:,range_max); %1024个chrip的某列range bin的相位
            vital_data = angle_fft_last;
            range_index = range_max;
            figure
            plot((1:M/fs),range_index*rangeRes*ones(1,M/fs));
            xlabel('时间(s)');
            ylabel('距离(m)');
            title('人体所在位置');
        case 2
        %% 2.取每帧的最大距离门
            %求最大距离门随时间的变化
            range_profile_abs=abs(data_out);
            range_max_index=ones(M,1);
            range_second_index=ones(M,1);
            for i=1:M
                for j=1:N
                    if((j*rangeRes)<2.5 &&(j*rangeRes)>0.3)
                       if (range_profile_abs(j,i,1)>range_profile_abs(range_max_index(i),i,1))
                           range_second_index(i)=range_max_index(j); 
                           range_max_index(i)=j;
                       elseif (range_profile_abs(j,i,1)>range_profile_abs(range_second_index(i),i,1))
                           range_second_index(i)=j;
                        end
                    end
                end
            end        
            % 对距离门突变进行校正
%             for i=5:M-4
%                 if(abs(range_max_index(i)-mean(range_max_index(i-4:i+4)))>2)
%                     range_max_index(i)=range_max_index(i-1);
%                 end
%             end
            for i=1:M
                angle_fft_last(i)=angle_fft(i,range_max_index(i));
            end
            vital_data = angle_fft_last;
            range_index = range_max_index;
            figure
            plot((1/fs:1/fs:M/fs),range_index*rangeRes);
            xlabel('时间(s)');
            ylabel('距离(m)');
            title('人体所在位置');
        case 3
        %% 3.根据一段时间内的距离门单元众数选择
            %求最大距离门随时间的变化
            range_profile_abs=abs(data_out);
            range_max_index=ones(M,1);
            range_second_index=ones(M,1);
            for i=1:M
                for j=1:N
                    if((j*rangeRes)<2.5 &&(j*rangeRes)>0.3)
                       if (range_profile_abs(j,i,1)>range_profile_abs(range_max_index(i),i,1))
                           range_second_index(i)=range_max_index(j); 
                           range_max_index(i)=j;
                       elseif (range_profile_abs(j,i,1)>range_profile_abs(range_second_index(i),i,1))
                           range_second_index(i)=j;
                        end
                    end
                end
            end        
            % 对距离门突变进行校正
%             for i=5:M-4
%                 if(abs(range_max_index(i)-mean(range_max_index(i-4:i+4)))>2)
%                     range_max_index(i)=range_max_index(i-1);
%                 end
%             end
            time_divide=fs*t; 
            for i=1:time_divide:M  
                angle_fft_last_index(floor(i/time_divide)+1)=mode(range_max_index(i:i+time_divide-1));
            end
            for i=1:M-1
                angle_fft_last(i)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)+1));
            end
            angle_fft_last(M)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)));
            vital_data = angle_fft_last;
            range_index = angle_fft_last_index;
            figure;
            plot((1:t:M/fs),range_index*rangeRes);
            xlabel('时间(s)');
            ylabel('距离(m)');
            title('人体所在位置');
        case 4
        %% 4.根据一段时间内的非相干累积值选择
            time_divide = fs*t;
            range_index2 = zeros(1,M/time_divide);
            range_profile_abs=abs(data_out);
                for i = 1:time_divide:M 
                    range_profile_last2=zeros(N,1);
                    range_max2=1;
                    for j = 1:N
                       if((j*rangeRes)<2.5 &&(j*rangeRes)>0.3) % 限定了检测距离为0.5-2.5m
                            for k= i:i+time_divide-1
                                range_profile_last2(j) = range_profile_last2(j) + power(range_profile_abs(j,k),2);
                            end
                            if (range_profile_last2(j)>range_profile_last2(range_max2))
                                range_max2=j;
                            end
                        end
                    end 
                    range_index2(floor(i/time_divide)+1) = range_max2;
                    angle_fft_last(i:i+time_divide-1) = angle_fft(i:i+time_divide-1,range_max2);
                end
                vital_data = angle_fft_last;
                range_index = range_index2;
                figure;
                plot((1:t:M/fs),range_index*rangeRes);
                xlabel('时间(s)');
                ylabel('距离(m)');
                title('人体所在位置');
%         case 5
        % 5.根据信噪比自适应选择距离门.
            % time_list=[0.5,1,2,4,6,10,15,20,30];
            % SN = zeros(1,length(time_list));
            % rms_value = zeros(1,length(time_list));
            % % 每过一段时间检查是否发生位移
            % for t=1:length(time_list)
            %     time_divide=fs*time_list(t); 
            %     for i=1:time_divide:numChirps  
            %         angle_fft_last_index(floor(i/time_divide)+1)=mode(range_max_index(i:i+time_divide-1));
            %     end
            %     for i=1:numChirps-1
            %         angle_fft_last(i)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)+1));
            %     end
            %     angle_fft_last(numChirps)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)));
            %     SN(t) = snr(angle_fft_last,50,3);
            % end
            % figure
            % X = time_list;
            % plot(X,SN);hold on;
            % [max_value,time_best] = max(SN);
            % % 根据信噪比选出最好的分段时间
            %     time_divide=fs*time_list(time_best); 
            %     for i=1:time_divide:numChirps  
            %         angle_fft_last_index(floor(i/time_divide)+1)=mode(range_max_index(i:i+time_divide-1));
            %     end
            %     for i=1:numChirps-1
            %         angle_fft_last(i)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)+1));
            %     end
            %     angle_fft_last(numChirps)=angle_fft(i,angle_fft_last_index(floor(i/time_divide)));
    end
end