function [ output_grid_tmp ] = restrict_op( tmp_grid )
%RESTRICT restrict operator for the multigrid method
%%  
grid_size = size(tmp_grid,1)+1;



output_grid_tmp = tmp_grid(1:2:grid_size-3,1:2:grid_size-3).*(1/16);
output_grid_tmp = output_grid_tmp + tmp_grid(1:2:grid_size-3,2:2:grid_size-2).*(1/8);
output_grid_tmp = output_grid_tmp + tmp_grid(1:2:grid_size-3,3:2:grid_size-1).*(1/16);

output_grid_tmp = output_grid_tmp + tmp_grid(2:2:grid_size-2,1:2:grid_size-3).*(1/8);
output_grid_tmp = output_grid_tmp + tmp_grid(2:2:grid_size-2,2:2:grid_size-2).*(1/4);
output_grid_tmp = output_grid_tmp + tmp_grid(2:2:grid_size-2,3:2:grid_size-1).*(1/8);

output_grid_tmp = output_grid_tmp + tmp_grid(3:2:grid_size-1,1:2:grid_size-3).*(1/16);
output_grid_tmp = output_grid_tmp + tmp_grid(3:2:grid_size-1,2:2:grid_size-2).*(1/8);
output_grid_tmp = output_grid_tmp + tmp_grid(3:2:grid_size-1,3:2:grid_size-1).*(1/16);



end

