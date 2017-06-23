function [ int ] = int_2d( func )
%ITER 2d function

func_1d = @(x)(int_1d(@(y)(func(x,y))));

int = int_1d(func_1d);

end

