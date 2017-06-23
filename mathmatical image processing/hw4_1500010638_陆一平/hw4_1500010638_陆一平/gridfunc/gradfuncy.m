function[A]=gradfuncy(m,n)
% approximation of nabla_x operator
% full(gradfunc2(3,4))
lambda = sqrt(2)-1;
A = gradfuncy2(m,n).*lambda+gradfuncy1(m,n).*(1-lambda);
end


function [A]=gradfuncy2(m,n)
A = speye(m).*(0.5);
B = sparse(n,n);


for i = 1:n-1
    B(i+1,i)=-1;
    B(i,i+1)=1;
end

A=kron(B,A);

end

function [ A ] = gradfuncy1( m,n )
%GRADFUNCY :the operator which approximates nabla _y

A = sparse(m,m);
B = sparse(n,n);

for i = 1:m-1
    A(i,i+1)=0.25;
    A(i+1,i)=0.25;
end

for i = 1:n-1
    B(i+1,i)=-1;
    B(i,i+1)=1;
end


A=kron(B,A);
end


function flag = is_safe(x,y,m,n)
    flag = true;
    if x<=0 || x>m
        flag = false;
    end
    if y<=0|| y>n
        flag = false;
    end
end

function pos = get_pos(x,y,m,n)
    pos = x + (y-1)*m;
end
