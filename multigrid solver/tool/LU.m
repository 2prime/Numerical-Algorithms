function [L,U] = LU( A )
%% LU factorization
n=size(A,1);
grid_size = sqrt(n);
for k=1:n-1
    A(k+1:min(k+grid_size,n),k)=A(k+1:min(k+grid_size,n),k)/A(k,k);
    A(k+1:min(k+grid_size,n),k+1:n)=A(k+1:min(k+grid_size,n),k+1:n)-A(k+1:min(k+grid_size,n),k)*A(k,k+1:n);
end
L=speye(n,n)+tril(A,-1);
U=triu(A);

end
