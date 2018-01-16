function [ int ] = int_1d( fun )
%INT_1D 1d function integral

% a,b������������
a = 0;
b = 1;

%n gauss int
n=30;

syms x

p=sym2poly(diff((x^2-1)^(n+1),n+1))/(2^n*factorial(n));
tk=roots(p); % ����ڵ�
% �������ϵ��
Ak=zeros(n+1,1);
for i=1:n+1
    xkt=tk;
    xkt(i)=[];
    pn=poly(xkt);
    fp=@(x)polyval(pn,x)/polyval(pn,tk(i));
    Ak(i)=integral(fp,-1,1); % ���ϵ��
end
% ���ֱ�����������[a,b]�任��[-1,1]
xk=(b-a)/2*tk+(b+a)/2;
% �����������֮����ֺ�����ֵ
fx=zeros(n+1,1);
for i=1:n+1
    fx(i)=fun(xk(i))*(b-a)/2;
end
% �������ֵ
int=sum(Ak.*fx);

end

