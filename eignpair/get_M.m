function [ M ] = get_M( grid_size )
%GET_M get matrix M

diag_block = 1/12*sparse(diag(ones(1,grid_size-2),1))+1/12*sparse(diag(ones(1,grid_size-2),-1))+speye(grid_size-1).*1/2;
block = 1/12*speye(grid_size-1)+1/12*sparse(diag(ones(1,grid_size-2),1));


kron_matrix = sparse(diag(ones(1,grid_size-2),-1));

M = kron(speye(grid_size-1),diag_block)+kron(kron_matrix',block)+kron(kron_matrix,block');


end

