function [x, out] = l1_cvx_mosek(x0, A, b, mu, opts)
cvx_begin
%   default mosek solver in cvx package calls for a professional license
%   cvx_solver mosek
    cvx_solver mosek_2
    x = x0;
    variable x(size(x, 1))
    minimize (0.5 * (A * x - b)' * (A * x - b) + mu * sum(abs(x)));
cvx_end