function [ result_grid ] = V_cyc_heat( grid_data,time_step )
%V_CYC_HEAT V cycle multi_grid
%%
%pay_attenion to the gride_data
%grid_data should not have the zero padding.
%%


%%
grid_size = size(grid_data,1)+1;

if grid_size > 8
    grid = im_to_vec(grid_data);
    Iter_opt = Iter_tool;
    Iter_opt.iter_time=3;
    A=implicit_matrix(grid_size,time_step );
    result_grid = g_s(A,grid,Iter_opt);
    res = grid-A*result_grid;
    
    %multi-grid solver for A*x=res
    res_new = restrict_op(vec_to_im(res));
    size(res_new)
    result_res = V_cyc_heat(res_new,time_step);
    result_res = im_to_vec(lifting_op(result_res));
    
    result_grid = result_grid + g_s_x0(A,res,result_res,Iter_opt);
    result_grid = vec_to_im(result_grid);
else
    grid = im_to_vec(grid_data);
    A = implicit_matrix(grid_size,time_step );
    result_grid = cgs(A,grid);
    result_grid = vec_to_im(result_grid);

end

end

