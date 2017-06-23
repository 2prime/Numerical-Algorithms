function [ output ] = explicit_heat( init,n,t )
%EXPLICIT_HEAT �ȷ�����ʾ��ʽ
%   ���룺init ��ʼֵ
%        n ʱ��������Ŀ
%        t  output_time

%%
time_step = t/n;

output = init;

for i = 1:n
    lap = lap_z(output);
    output = output + time_step * lap;
end


end

