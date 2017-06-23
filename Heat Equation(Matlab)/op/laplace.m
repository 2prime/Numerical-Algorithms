function [ output_func ] = laplace( func )
%LAPLACE get the laplace of func
%   laplace = u_xx+u_yy
%   five point finite difference scheme
%%
lap = [0,1,0;1,-4,1;0,1,0];

grid_cof = (size(func,1)-1)^2;

output_func = conv2(func,lap*grid_cof,'valid');


end

