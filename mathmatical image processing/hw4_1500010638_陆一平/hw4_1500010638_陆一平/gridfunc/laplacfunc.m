function [ A ] = laplacfunc( m,n )
%LAPLACFUNC : approximate for laplace operator
% tmp = sparse(ones(m,m).*(1/3));
% B = spdiags(ones(n,1),-1,n,n)+spdiags(ones(n,1),1,n,n);
% 
% A = kron(B,tmp) + kron(speye(n),tmp-eye(m));

tmp1 = kron(spdiags(ones(m, 1), -1, m, m), 1/3) + kron(spdiags(ones(m, 1), 1, m, m), 1/3) + kron(spdiags(ones(m, 1), 0, m, m), 1/3);
tmp1(1, m) = 1/3;tmp1(m, 1) = 1/3;

tmp2 = kron(spdiags(ones(m, 1), 1, m, m), 1/3) + kron(spdiags(ones(m, 1), -1, m, m), 1/3) + kron(spdiags(ones(m, 1), 0, m, m), -8/3);
tmp2(1, m) = 1/3;tmp2(m, 1) = 1/3;

A = kron(spdiags(ones(n, 1), 1, n, n), tmp1) + kron(spdiags(ones(n, 1), -1, n, n), tmp1) + kron(spdiags(ones(n, 1), 0, n, n), tmp2);
 
A(m * (n - 1) + 1: m * n, 1 : m) = tmp1;
A(1 : m, m * (n - 1) + 1: m * n) = tmp1;

end

