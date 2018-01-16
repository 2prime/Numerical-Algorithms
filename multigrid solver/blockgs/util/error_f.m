function e_f=error_f(x,fx)
%%
% Error function

N=size(x,1)-1;
xx=interpolate(x);
xx0=xx(2:2:2*N,2:2:2*N);
row=0:1/(2*N):1;
row=row(1:2*N+1);
Col=ones(2*N+1,1)*row;
Col=Col(2:2:2*N,2:2:2*N);
FXX=arrayfun(fx,Col',Col);
e_f=err(xx0-FXX);

end
