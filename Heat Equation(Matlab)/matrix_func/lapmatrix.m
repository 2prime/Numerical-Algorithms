function [ lap_matrix ] = lapmatrix( grid_size )
%LAPMATRIX get the matrix of the laplace operator

diag_block = sparse(diag(ones(1,grid_size-2),1))+sparse(diag(ones(1,grid_size-2),-1))+speye(grid_size-1).*-4;
block = speye(grid_size-1);
%block(1,grid_size-1)=1;

kron_matrix = sparse(diag(ones(1,grid_size-2),-1));

lap_matrix=kron(speye(grid_size-1),diag_block)+kron(kron_matrix,block)+kron(kron_matrix',block');

lap_matrix = lap_matrix.*(grid_size^2);

end

