function [x,err_tol,iter]=cg_aniso(rhs,eps,grid_size)

x=ones((grid_size-1)^2,1);
%x=zeros((grid_size-1)^2,1);
A = anisolap(grid_size,eps);
b = rhs;
maxit = 200;


r = b - A*x;
err = norm(r);
end_tol = 1e-8*err;
err_tol = zeros(1,maxit);
iter = 0;


while (iter < maxit) && (err>end_tol)
    iter = iter + 1;
    z = r;
    if iter == 1
        p = z;
    else

        beta = (r'*z)/(old_r'*old_z);
        p = z + beta*p;
    end
    alpha = (r'*z)/(p'*(A*p));
    x = x + alpha * p;
    old_r = r;
    old_z = z;
    r = r-alpha*(A*p);
    err = norm(r);
    err_tol(iter)=err;
end




end
