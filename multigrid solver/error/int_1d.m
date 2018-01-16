function [ int ] = int_1d( fun )
%INT_1D 1d function integral

% a,b：积分上下限
a = 0;
b = 1;

%n gauss int
n=30;

syms x

p=sym2poly(diff((x^2-1)^(n+1),n+1))/(2^n*factorial(n));
tk=roots(p); % 求积节点
% 计算求积系数
Ak=zeros(n+1,1);
for i=1:n+1
    xkt=tk;
    xkt(i)=[];
    pn=poly(xkt);
    fp=@(x)polyval(pn,x)/polyval(pn,tk(i));
    Ak(i)=integral(fp,-1,1); % 求积系数
end
% 积分变量代换，将[a,b]变换到[-1,1]
xk=(b-a)/2*tk+(b+a)/2;
% 计算变量代换之后积分函数的值
fx=zeros(n+1,1);
for i=1:n+1
    fx(i)=fun(xk(i))*(b-a)/2;
end
% 计算积分值
int=sum(Ak.*fx);

end

