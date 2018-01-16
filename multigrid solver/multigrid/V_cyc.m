function [ result_grid ] = V_cyc( grid_size,eps,rhs )
%V_CYC_HEAT V cycle multi_grid
%%
%pay_attenion to the gride_data
%grid_data should not have the zero padding.
%%


%%

if grid_size > 8
    Iter_opt = Iter_tool;
    Iter_opt.iter_time=3;
    A=anisolap(grid_size,eps);
    result_grid = g_s(A,rhs,Iter_opt);
    res = rhs-A*result_grid;
    
    %multi-grid solver for A*x=res
    res_new = restrict_op(vec_to_im(res));
    size(res_new)
    result_res = V_cyc(size(res_new,1)+1,eps,im_to_vec(res_new));
    result_res = im_to_vec(lifting_op(result_res));
    
    result_grid = result_grid + g_s_x0(A,res,result_res,Iter_opt);
    result_grid = vec_to_im(result_grid);
else
    %small grid: solve directly
    A = anisolap(grid_size, eps);
    result_grid = cgs(A,rhs);
    result_grid = vec_to_im(result_grid);
    size(result_grid)

end

end

