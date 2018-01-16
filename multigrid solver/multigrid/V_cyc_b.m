function [ result_grid ] = V_cyc( grid_size,eps,rhs )
%V_CYC_HEAT V cycle multi_grid
%%
%grid_data should not have the zero padding.
%%


%%

if grid_size > 8
    Iter_opt = Iter_tool;
    Iter_opt.iter_time=3;
    A=anisolap(grid_size,eps);
    result_grid = GSB(rhs,eps);
    res = rhs-A*result_grid;
    
    %multi-grid solver for A*x=res
    res_new = restrict_op(vec_to_im(res));
    size(res_new)
    result_res = V_cyc(size(res_new,1)+1,eps,im_to_vec(res_new));
    result_res = im_to_vec(lifting_op(result_res));
    
    result_grid = result_grid + GSB_x0(result_res,res,eps);
    result_grid = vec_to_im(result_grid);
else
    A = anisolap(grid_size, eps);
    result_grid = cgs(A,rhs);
    result_grid = vec_to_im(result_grid);
    size(result_grid)

end

end

