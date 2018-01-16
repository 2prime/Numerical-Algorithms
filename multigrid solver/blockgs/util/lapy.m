function r = lapy(u,h,e);
%% u_yy matrix

[nx,ny] = size( u );
r = zeros(nx,ny);
h2 = h*h;
r(2:nx-1,2:ny-1) =(e*u(1:nx-2,2:ny-1)+e*u(3:nx,2:ny-1)...
    +u(2:nx-1,1:ny-2)+u(2:nx-1,3:ny)-(2+2*e)*u(2:nx-1,2:ny-1))/h2;
end
