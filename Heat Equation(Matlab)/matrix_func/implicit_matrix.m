function [ A ] = implicit_matrix( grid_size,time_step )
%IMPLICIT_MATRIX 此处显示有关此函数的摘要
%   此处显示详细说明

lap_matrix = lapmatrix(grid_size).*time_step;
A = speye((grid_size-1)^2)-lap_matrix;

end

