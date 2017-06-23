function phi=Reinitial2D(phi,steps)
%%%%reinitialization
[M,N]=size(phi); h=1/(max(N,M)-1);eps=1e-9;
CFL=0.5;
dt=h*CFL; dx=1/(M-1); dy=1/(N-1);
for k=1:steps
    xp = (phi([2:end,end],:)-phi)/dx;
    xn = (phi-phi([1,1:end-1],:))/dx;
    yp = (phi(:,[2:end,end])-phi)/dy;
    yn = (phi-phi(:,[1,1:end-1]))/dy;
    phi_p = phi>=0;
    phi_n = 1-phi_p;
    Godnov_p = sqrt(max((max(xn,0)).^2,(min(xp,0)).^2)+max((max(yn,0)).^2,(min(yp,0)).^2));
    Godnov_n = sqrt(max((min(xn,0)).^2,(max(xp,0)).^2)+max((min(yn,0)).^2,(max(yp,0)).^2));
    phi = phi - dt*phi_p.*(Godnov_p-1).*phi./sqrt(phi.^2+Godnov_p.^2*dx*dy+eps) ...
              - dt*phi_n.*(Godnov_n-1).*phi./sqrt(phi.^2+Godnov_n.^2*dx*dy+eps);
end