function [ error ] = get_error( grid1,t )
%GET_ERROR 此处显示有关此函数的摘要
%   此处显示详细说明

error_func = @(x,y)(get_value(grid1,x,y)-exp(-2*t*(pi^2))*sin(pi*x)*sin(pi*y));

error_func_l2 = @(x,y)(error_func(x,y)^2);

error = int_2d(error_func_l2);

end

