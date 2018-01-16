function u = gsone(f,u,h,e)

[nx,ny] = size( f );
h2 = h*h;
r2 = 0;
T=1/(2+2*e);
for j=2:ny-1
  for i=2:nx-1
    u(i,j) = (e*u(i-1,j) + e*u(i+1,j) + u(i,j-1) + u(i,j+1) + h2*f(i,j))*T;
  end
end
end
