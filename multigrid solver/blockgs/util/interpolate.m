function u = interpolate( uc )
%% lifiting operator
[nx,ny] = size( uc );
nxf = 2*(nx-1)+1;
nyf = 2*(nx-1)+1;
u = zeros(nxf,nyf);

for j=1:ny
  for i=1:nx
    u(2*(i-1)+1,2*(j-1)+1) = uc(i,j);
  end
end

for j=1:2:nyf
  for i=2:2:nxf-1
    u(i,j) = 0.5*(u(i-1,j) + u(i+1,j));
  end
end

for j=2:2:nyf-1
  for i=1:nxf
    u(i,j) = 0.5*(u(i,j-1) + u(i,j+1));
  end
end
