function [c g]=ChebyshevApprox(f,n)
%%%% Assuming f : [0 pi] -> R.
c=zeros(1,n);
a=pi/2;
for k=1:n
    Integrand = @(x)(cos((k-1)*x).*f(a*(cos(x)+1)));
    c(k)=2/pi*quad(Integrand,0,pi);
end
g=@(x)(ChebyshevPoly(x,c,a));

function p=ChebyshevPoly(x,c,a)
n=length(c);
% p=1/2*c(1);
% for k=2:n
%     p=p+c(k)*cos((k-1)*acos((x-a)/a));
% end
T0=1;T1=(x-a)/a;
p=1/2*c(1)*T0+c(2)*T1;
for k=3:n
    Tk=2/a*(x-a).*T1-T0;
    T0=T1;T1=Tk;
    p=p+c(k)*Tk;
end

