function e_g=error_g(x,g1,g2)
%%
% Error Of g
N=size(x,1)-1;
h=1/(2*N);
xx=interpolate(x);
x_grad=zeros(N,N,2);
for i=1:N
    for j=1:N
        x_grad(i,j,:)=[(xx(2*i+1,2*j)-xx(2*i-1,2*j))/(2*h),(xx(2*i,2*j+1)-xx(2*i,2*j-1))/(2*h)];
    end
end
row=0:1/(2*N):1;
row=row(1:2*N+1);
Col=ones(2*N+1,1)*row;
Col=Col(2:2:2*N,2:2:2*N);
GX(:,:,1)=arrayfun(g1,Col',Col);
GX(:,:,2)=arrayfun(g2,Col',Col);
e_g=err(max(GX-x_grad));
end
