function [x,info]=GS(x0,f,e)

x=x0;
N=size(x,1)-1;
h=1/N;
loss0=1;
loss=loss0;
k=0;
while(1)
k=k+2;
x=gsone(f,x,h,e);
x=gsinv(f,x,h,e);
    loss=err(f+lapy(x,h,e));
    if loss/loss0<10^-8
        break;
    end
end
info={};
info.iter=k;
end
