function [x,info]=GSB(x0,f,e)
x=x0;
N=size(x,1)-1;
h=1/N;
loss0=1;
loss=loss0;
k=0;
n=N-1;
T=2*eye(n);
I=eye(n);
T(1:n-1,2:n)=T(1:n-1,2:n)-eye(n-1);
T(2:n,1:n-1)=T(2:n,1:n-1)-eye(n-1);
P=(2*e*I+T)^-1;
while(1)
k=k+2;
x=gsbone(f,x,h,e,P); %正向更新一遍块
x=gsbinv(f,x,h,e,P); %反向更新一边块
loss=err(f+lapy(x,h,e));
if loss/loss0<10^-8
    break;
end
end
info={};
info.iter=k;
end
