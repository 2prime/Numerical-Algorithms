function x3=test(N,e)
%%
% N grid size

func=@(x,y)((1+e)*pi^2*sin(pi*x)*sin(pi*y)); % func on the right
fx=@(x,y)(sin(pi*x)*sin(pi*y)); % Truth Solution
g1=@(x,y)(pi*cos(pi*x)*sin(pi*y));% u_xx
g2=@(x,y)(pi*sin(pi*x)*cos(pi*y));% u_yy
h=1/N;

%% BGS Method
row=0:h:1;
row=row(1:N+1);
Col=ones(N+1,1)*row;

F=arrayfun(func,Col',Col);
FX=arrayfun(fx,Col',Col);

b=F;
x0=zeros(N+1,N+1);
x0(2:N,2:N)=rand(N-1,N-1);

fprintf('\nN=%d,e=%d\n',N,e)

%%
tic;
[x1,info1]=GS(x0,b,e);
t1=toc;
fprintf('G-S Method Iter:%6d       time:%f   e_f:%2.2d   e_g:%2.2d\n',info1.iter,t1,error_f(x1,fx),error_g(x1,g1,g2))

%%
tic;
[x2,info2]=GSB(x0,b,e);
t2=toc;
fprintf('Block G-S Iter:%6d       time:%f   e_f:%2.2d   e_g:%2.2d\n',info2.iter,t2,error_f(x2,fx),error_g(x2,g1,g2))
end
