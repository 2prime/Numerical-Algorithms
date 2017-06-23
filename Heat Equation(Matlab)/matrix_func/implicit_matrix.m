function [ A ] = implicit_matrix( grid_size,time_step )
%IMPLICIT_MATRIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

lap_matrix = lapmatrix(grid_size).*time_step;
A = speye((grid_size-1)^2)-lap_matrix;

end

