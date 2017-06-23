function [ answer ] = test_answer( t,grid_size )
%TEST_ANSWER get the ture answer at time t
%   true answer of the equation: u_t = Laplace u

%%
x = 0:1/grid_size:1;
y = 0:1/grid_size:1;

[X,Y] = meshgrid(x,y);


answer = exp(-2*t*(pi^2))*sin(pi*X).*sin(pi*Y);

end

