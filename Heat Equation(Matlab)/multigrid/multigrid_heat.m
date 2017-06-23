function [ result_grid ] = multigrid_heat( grid_data,time_step,v_num)
%MULTIGRID 此处显示有关此函数的摘要
%   此处显示详细说明

grid_size = size(grid_data,1)+1;
A=implicit_matrix(grid_size,time_step );
func_tmp = @(x)(vec_to_im(A*(im_to_vec(x))));

result_grid = V_cyc_heat(grid_data,time_step);
for i = 1:v_num-1
    residual_grid = grid_data - func_tmp(result_grid);
    result_grid = result_grid + V_cyc_heat(residual_grid,time_step);
end

end

