function [ output ] = explicit_heat( init,n,t )
%EXPLICIT_HEAT 热方程显示格式
%   输入：init 初始值
%        n 时间网格数目
%        t  output_time

%%
time_step = t/n;

output = init;

for i = 1:n
    lap = lap_z(output);
    output = output + time_step * lap;
end


end

