function [ output_grid ] = lifting_op( input_grid )
%LIFTING lifting operator for multigrid solver

%%
grid_size = size(input_grid,1)+1;
%add the zero padding
grid = zeros(grid_size+1,grid_size+1);
grid(2:grid_size,2:grid_size)=input_grid;

output_grid = zeros(grid_size*2+1,grid_size*2+1);
output_grid(1:2:end,1:2:end) = grid;
output_grid(2:2:end,2:2:end) = (grid(1:grid_size,1:grid_size)+grid(1:grid_size,2:grid_size+1)+grid(2:grid_size+1,1:grid_size)+grid(2:grid_size+1,2:grid_size+1))./4;
output_grid(1:2:end,2:2:end) = (grid(:,1:grid_size)+grid(:,2:grid_size+1))./2;
output_grid(2:2:end,1:2:end) = (grid(1:grid_size,:)+grid(2:grid_size+1,:))./2;

output_grid = output_grid(2:grid_size*2,2:grid_size*2);

end

