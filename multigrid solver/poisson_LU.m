clear;clc;
addpath('op')
addpath('pdefunc')
addpath('plot_func')
addpath('matrix_func')
addpath('tool')
addpath('error')
addpath('multigrid')



%%
% parameter setting
grid_size = 32;
eps = 1e-2;



%%
% solution generating
x = 0:1/grid_size:1;
y = 0:1/grid_size:1;
[X,Y] = meshgrid(x,y);

answer = sin(pi*X).*sin(pi*Y); % imshow(answer);
size(answer(2:grid_size,2:grid_size)) % remove the zero padding
answer = im_to_vec(answer(2:grid_size,2:grid_size));
rhs = -(1+eps)*pi*pi*answer; % calculate the right hand side


lapm = anisolap(grid_size,eps);

tic;

[L,U] = LU(lapm);

%%%%%%%%%%%%
% LUx = rhs
%%%%%%%%%%%%

y = L\rhs;
x = U\y;
solu = zeros(grid_size+1);
solu(2:grid_size,2:grid_size) = vec_to_im(x);

imshow(vec_to_im(x));
fprintf('===================================\nInformation Print\n')
toc
fprintf('Vector Error:%g\n',norm(x-answer));
fprintf('Function Error:%g\n',get_error(solu, 0))