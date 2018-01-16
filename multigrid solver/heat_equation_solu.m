clc;
clear;
addpath('op')
addpath('pdefunc')
addpath('plot_func')
addpath('matrix_func')
addpath('tool')
addpath('error')
addpath('multigrid')


%%
%%NLA Func Demo

Iter_opt = Iter_tool;
Iter_opt.iter_time = 1000;
Iter_opt.iter_end = 1e-8;

msize = 3;
X = randn(msize);
A = X'*X;
b = randn(msize,1);

x = g_s(A,b,Iter_opt);

norm(A*x-b)

%%
%%Integral tool box Demo
func_test = @(x,y)(sin(pi*x)*sin(pi*y));
int_2d(func_test)
func_test = @(x)(sin(pi*x));
integral(func_test,0,1)^2


%% PDE Func Demo
grid_size = 128;
time = 0.2;
step = 10;
init = test_answer(0,grid_size);
tic;
tool=Scheme_tool;
tool.method='pcg';
output_cn = CN_heat(init,step,time,tool);

%%%implicit scheme
tic;
tool=Scheme_tool;
tool.method='pcg';
output_i = implicit_heat(init,step,time,tool);

fprintf('CN scheme step:%g,time:%g, grid error:%g\n',step,toc,norm(output_cn-test_answer(time,grid_size)));
fprintf('CN scheme eroor:%g\n',get_error(output_cn,time))
fprintf('implicit scheme 2norm,step:%g,time:%g,grid error:%g\n',step,toc,norm(output_i-test_answer(time,grid_size)));
fprintf('implicit scheme 2norm eroor:%g\n',get_error(output_i,time))

%%%explicit shcheme
step = 10;
tic;
output_e = explicit_heat(init,step,time);
fprintf('explicit scheme 2norm,step:%g,time:%g, error:%g\n',step,toc,norm(output_e-test_answer(time,grid_size)));
fprintf('implicit scheme 2norm eroor:%g\n',get_error(output_e,time))

%% demo vedio

grid_size = 128;
step_size = 2e-3;
answer_E = @(x,grid_size)explicit_heat(test_answer(0,grid_size),floor(x/step_size),x,tool);
answer_CN = @(x,grid_size)CN_heat(test_answer(0,grid_size),floor(x/step_size),x,tool);
error_CN = @(x,grid_size)abs(answer_CN(x,grid_size)-test_answer(x,grid_size));
tool.method='multigrid';
answer_i = @(x,grid_size)implicit_heat(test_answer(0,grid_size),floor(x/step_size),x,tool);



%show_surf(grid_size,error_CN,'Demo of the error(CN Scheme)','gridError',4e-5);
%bad_answer_E = @(x,grid_size)explicit_heat(test_answer(0,grid_size),100,x,tool);
%show_surf(grid_size,bad_answer_E,'Demo of heat equation(Explict Scheme Bad Answer)','Explict Scheme Bollow up',1);
%show_surf(grid_size,answer_E,'Demo of heat equation(Explict Scheme Answer)','Explict Scheme',1);
show_surf(grid_size,answer_CN,'Demo of heat equation(CN Scheme)','CN Scheme',1);
%show_surf(grid_size,answer_i,'Demo of heat equation(implicit Scheme)','implicit Scheme',1);
%show_surf(grid_size,@test_answer,'Demo of heat equation(True Answer)','True Answer',1);


