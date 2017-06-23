%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image toolbox
% MATLAB file
% 
% PeronaMalikExample.m
% Perona Malik anisotropic filer
%
% function B = PeronaMalik( A, b, tfinal, outfilename )
%
%
%   u' = div[ g(u) grad(u) ]
%   du/dn |_{d Omega} = 0
%   u(t=0) = u0
%
%   g(u)= 1/ (1+b |grad(u)|^2)
%
% input:  A: image file name
%         b: parameter in the edge detector
%         tfinal: artificial time
%         outfilename: output file name
% output: B: denoised image
% example: B = gaussian( A, 0.1, 10 );
%
% created:       07/28/2009
% author:        syleung@gmail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = PeronaMalik( A, b, tfinal )

fileout = 5;    % number of intermediate output

dtout = tfinal/fileout;
time=0;
[m,n,p]=size(A);
xi=2:m+1;
yi=2:n+1;
B=zeros(m+2,n+2);
B(xi,yi) = A(:,:,1);

while time < tfinal
    
    % du/dn=0 on the boundary
    
    B(1,yi)= B(2,yi);
    B(end,yi)= B(end-1,yi);
    B(xi,1)= B(xi,2);
    B(xi,end)= B(xi,end-1);
        
    % g(|grad u|^2)
    
    %C = gaussian(B,2,5);
    [gn,ge,gs,gw]=g(B,b);
    
    Bn = B(xi,yi+1);
    Be = B(xi+1,yi);
    Bs = B(xi,yi-1);
    Bw = B(xi-1,yi);
    Bc = B(xi,yi);
    
    % u^k+1 = u^k + dt * div[ g(u) grad(u) ]
    
    temp = gn.*Bn + ge.*Be + gs.*Bs + gw.*Bw ...
        - (gn+ge+gs+gw).*Bc;
    
    dt = 0.1;  
    time=time+dt;
    [time dt tfinal];
    
    B(xi,yi)=B(xi,yi)+dt*temp;
    
    
    
end

B=B(2:end-1,2:end-1);

return


function [gn,ge,gs,gw]=g(u,b)

% g(u)= 1/ (1+b |grad(u)|^2)

[m,n]=size(u);
xi=[2:m-1];
yi=[2:n-1];

uxp=u(xi+1,yi);
uxm=u(xi-1,yi);
uyp=u(xi,yi+1);
uym=u(xi,yi-1);
uc=u(xi,yi);

duxp=uxp-uc;
duxm=uc-uxm;
duyp=uyp-uc;
duym=uc-uym;

duxc=(duxp+duxm)/2;
duyc=(duyp+duym)/2;

gn=1./(1+b*(duxc.^2+duyp.^2));
ge=1./(1+b*(duxp.^2+duyc.^2));
gs=1./(1+b*(duxc.^2+duym.^2));
gw=1./(1+b*(duxm.^2+duyc.^2));

return