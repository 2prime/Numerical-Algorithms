function [ U ] = eignpair( A,M )
%EGINPAIR Ax=lambda Mx
[w,h]=size(A);
U = ones(w,1);
U = U/sqrt((M*U)'*U);
rou_0 = ((A*U)'*U)/((M*U)'*U);
rou = rou_0;


for k = 1:300
    Z = orthbase(A,M,rou,U,10);
    A_0 = Z'*(A-rou*M)*Z;
    M_0 = Z'*M*Z;
    [V,D]=eigs(A_0,M_0,1,'sm');
    rou = rou+D;
    U = Z*V;
    U = U/sqrt((M*U)'*U);
    if norm(A*U-rou*M*U)<1e-3
        break
    end
end

fprintf('iter time:%g\n',k);
fprintf('Error:%g\n',norm(A*U-rou*M*U))

end

