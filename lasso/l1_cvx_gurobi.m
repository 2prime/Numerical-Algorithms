function [x, out] = l1_cvx_gurobi(x0, A, b, mu, opts)
cvx_begin quiet
%   default gurobi solver in cvx package calls for a professional license
%   cvx_solver gurobi
    cvx_solver gurobi_2
    x = x0;
    variable x(size(x))
    minimize (0.5 * (A * x - b)' * (A * x - b) + mu * sum(abs(x)));
cvx_end