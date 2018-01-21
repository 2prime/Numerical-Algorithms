function [ lap_matrix ] = anisolap( grid_size, eps)
%ANISOLAP
%%sparse matrix of anisotropic laplace

diag_block = eps*sparse(diag(ones(1,grid_size-2),1))+eps*sparse(diag(ones(1,grid_size-2),-1))+speye(grid_size-1).*-2*(1+eps);
block = speye(grid_size-1);
%block(1,grid_size-1)=1;

kron_matrix = sparse(diag(ones(1,grid_size-2),-1));

lap_matrix=kron(speye(grid_size-1),diag_block)+kron(kron_matrix,block)+kron(kron_matrix',block');

lap_matrix = lap_matrix.*(grid_size^2);


end

