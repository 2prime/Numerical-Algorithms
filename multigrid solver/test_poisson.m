%%%%%
% Code for Report
% Numerical Scheme For Anisotropic Poisson Equation With First Boundary Condition On Square Domain
% Yiping Lu Peking Univ
% @2prime 2018/1
%%%%%

clc;clear;
addpath('op')
addpath('pdefunc')
addpath('plot_func')
addpath('matrix_func')
addpath('tool')
addpath('error')
addpath('multigrid')

%%
% parameter setting
grid_size = 128;
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

Iter_opt = Iter_tool;
Iter_opt.iter_time=100;

%% get the solution of the problem
answer_V = V_cyc(grid_size,eps,rhs);
subplot(1,4,1);
imshow(answer_V);
title('solution via multigrid');
subplot(1,4,2);
tic;
[pcg_solu,err_pcg,pcg_it] = pcg_aniso(rhs,eps,grid_size);
toc
imshow(vec_to_im(pcg_solu));
title('solution via pcg');
subplot(1,4,3);
[cg_solu,err_cg,cg_it] = cg_aniso(rhs,eps,grid_size);
imshow(vec_to_im(cg_solu));
title('solution via cg');
subplot(1,4,4);
imshow(vec_to_im(answer));
title('true solution');


% answer_V = g_s(lapm,rhs,Iter_opt);
% imshow(vec_to_im(answer_V));

figure;
semilogy(err_pcg);
hold on;
semilogy(err_cg);
legend('pcg','cg');
ylabel('Error');
xlabel('Iteration Time');

solu = zeros(grid_size+1);
solu(2:grid_size,2:grid_size) = vec_to_im(pcg_solu);
fprintf('=====================================\n Information printing:\n');
fprintf('grid size:%g, eps:%g\n',grid_size,eps);
fprintf('pcg iteration times:%g\n',pcg_it);
fprintf('vector l2 norm error:%g\n',norm(pcg_solu-answer));
fprintf('function l2 norm error:%g\n',get_error(solu, 0));



