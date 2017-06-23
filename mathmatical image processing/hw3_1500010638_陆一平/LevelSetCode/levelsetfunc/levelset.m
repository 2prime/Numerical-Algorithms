function [v,g]=levelset(phi,steps,renstep)
%% levelset method
%authorrized by @2pimre
%2017.3
[m,n,p]=size(phi);
gamma= 2;
alpha = 0.5;
beta= 4.5;
step_size = 0.05;

%% preprocessing for picture
%calculate g
xi=2:m+1;
yi=2:n+1;
B=zeros(m+2,n+2);%image
B(xi,yi) = phi(:,:,1);
B(1,yi)= B(2,yi);
B(end,yi)= B(end-1,yi);
B(xi,1)= B(xi,2);
B(xi,end)= B(xi,end-1);

%calculate g(|\nabla I|)
%where g(s) = 1/(1+s^2)
g=cal_g(B);
[dg_x,dg_y,~,~,~] =grad(g);



%% begin levelset interate

v=zeros(m,n);%level set function
for i = 1:m
    for j = 1:n
        v(i,j) = -(min([i j m-i n-j])-10);
    end
end

iter_time = 0;
for iter = 1:steps
    [vx,vy,vxx,vyy,vxy]=grad(v);
    [vx1,vx2,vy1,vy2]=two_grad(v);
    
    meancurv = ((vx.^2).*vyy+(vy.^2).*vxx-2*vx.*vy.*vxy)./(max((vx.^2+vy.^2),1e-10)).*g;
    
    evolution = ((max(vx1,0).^2+min(vx2,0).^2+max(vy1,0).^2+min(vy2,0).^2).^(1/2)).*g;
    
    advection = (max(dg_x,0).*vx1+min(dg_x,0).*vx2+max(dg_y,0).*vy1+min(dg_y,0).*vy2);
    
    v = v+ step_size*(gamma*meancurv+alpha*evolution+beta*advection);
    
    iter_time = iter_time+1;
    if iter_time == renstep
        v = Reinitial2D(v,20);
        iter_time = 0;
    end
    
end

v = Reinitial2D(v,20);



end

function [g]=cal_g(u)
%% g(u)= 1/ (1+ |grad(u)|^2)

[m,n]=size(u);
xi=[2:m-1];
yi=[2:n-1];

uxp=u(xi+1,yi);
uxm=u(xi-1,yi);
uyp=u(xi,yi+1);
uym=u(xi,yi-1);


%calculate grad
duxc=(uxp-uxm)/2;
duyc=(uyp-uym)/2;

g=1./(1+0.5*(duxc.^2+duyc.^2));


end

function [gx1,gx2,gy1,gy2]=two_grad(func)
%% calculate grad
    [m,n]=size(func);
    B=zeros(m+2,n+2);
    xi=2:m+1;
    yi=2:n+1;
    B(xi,yi) = func;
    B(1,yi)= B(2,yi);
    B(end,yi)= B(end-1,yi);
    B(xi,1)= B(xi,2);
    B(xi,end)= B(xi,end-1);
    B(1,1)=B(2,2);
    B(1,end)=B(1,end-1);
    B(end,end)=B(end-1,end-1);
    B(end,1)=B(end-1,1);
    
    uxp=B(xi+1,yi);
    uxm=B(xi-1,yi);
    uyp=B(xi,yi+1);
    uym=B(xi,yi-1);
    uc=B(xi,yi);
    
    %calculate grad
    gx1=uxp-uc;gx2=uc-uxm;
    gy1=uyp-uc;gy2=uc-uym;
    
    
    
end


function [gx,gy,gxx,gyy,gxy]=grad(func)
%% calculate grad
    [m,n]=size(func);
    B=zeros(m+2,n+2);
    xi=2:m+1;
    yi=2:n+1;
    B(xi,yi) = func;
    B(1,yi)= B(2,yi);
    B(end,yi)= B(end-1,yi);
    B(xi,1)= B(xi,2);
    B(xi,end)= B(xi,end-1);
    B(1,1)=B(1,2);
    B(1,end)=B(1,end-1);
    B(end,end)=B(end-1,end-1);
    B(end,1)=B(end-1,1);
    
    uxp=B(xi+1,yi);
    uxm=B(xi-1,yi);
    uyp=B(xi,yi+1);
    uym=B(xi,yi-1);
    uc=B(xi,yi);
    uxpyp=B(xi+1,yi+1);
    uxpym=B(xi+1,yi-1);
    uxmyp=B(xi-1,yi+1);
    uxmym=B(xi-1,yi-1);
    %calculate grad
    gx=(uxp-uxm)/2;
    gy=(uyp-uym)/2;
    gxx=(uxp+uxm-2*uc);
    gyy=(uyp+uym-2*uc);
    gxy=(uxpyp+uxmym-uxpym-uxmyp)/4;
    
end





