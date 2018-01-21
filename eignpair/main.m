clc;clear;
grid_size = 64;
eps = 1e-5;

A = anisolap(grid_size,eps);
M = get_M(grid_size);

U = eignpair(A,M);

answer = zeros(grid_size+1,grid_size+1);
answer(2:grid_size,2:grid_size)=vec_to_im(U);

figure;
x = 0:1/grid_size:1;
y = 0:1/grid_size:1;

[X,Y] = meshgrid(x,y);

surf(X,Y,answer);
