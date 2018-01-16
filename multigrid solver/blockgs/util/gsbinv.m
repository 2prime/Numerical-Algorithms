function u = gsbinv(f,u,h,e,P)
% P = inv(2*e*I+T)
[nx,ny] = size( f );
h2 = h*h;
 n=nx-2;
for i=ny-1:2
    u(i,2:n+1) = (e*u(i-1,2:n+1) + e*u(i+1,2:n+1) + h2*f(i,2:n+1))*P;
end
end
