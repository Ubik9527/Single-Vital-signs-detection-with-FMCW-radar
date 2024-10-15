%% 在数组的某一个区间内寻找峰值
% 输入参数：
% x：数组
% x_start：数组起始位置
% x_end：数组结束位置
% 输出参数：
% index：最大值索引（从头开始的位置）
% max_value：最大值
function [index,max_value] = findmax(x,x_start,x_end)
    max_value = 0;
    for i = x_start:x_end
            if(x(i)>max_value)
                max_value = x(i);
                index=i;
            end
    end
end