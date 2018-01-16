function [ value ] = get_value( grid_func,x,y )
%GET_VALUE 此处显示有关此函数的摘要
%   此处显示详细说明

grid_size = size(grid_func,1)-1;
grid_rad = 1/grid_size;

n1 = x/grid_rad;
n2 = y/grid_rad;
x = n1-floor(n1);
y = n2-floor(n2);
n1 = floor(n1);
n2 = floor(n2);

value = inter_sqaure(grid_func(n1+1:n1+2,n2+1:n2+2),x,y);


end

