function [x, out] = l1_mosek(x0, A, b, mu, opts)
clear prob;
[r, res] = mosekopt('symbcon');
m = size(A, 1);
n = size(A, 2);
%t, t1, t2, x, y, s
prob.c = [0.5, zeros(1, m + n + 2), mu * ones(1, n)];
prob.a = zeros(m + 2 * n + 2, m + 2 * n + 3);
prob.a(1 : m, 4 : m + n + 3) = [A, -0.5 * eye(m)];
prob.a(m + 1, 1) = 1;
prob.a(m + 1, 2) = -1;
prob.a(m + 2, 1) = 1;
prob.a(m + 2, 3) = -1;
for i = 1 : n
    prob.a(m + 2 + i, 3 + i) = 1;
    prob.a(m + 2 + i, m + n + 3 + i) = 1;
    prob.a(m + n + 2 + i, 3 + i) = -1;
    prob.a(m + n + 2 + i, m + n + 3 + i) = 1;
end
prob.a = sparse(prob.a);
prob.blc = [b', -1, 1, zeros(1, 2 * n)]';
prob.buc = [b', -1, 1, inf * ones(1, 2 * n)]';
prob.blx = [-inf * ones(n + m + 3, 1); zeros(n, 1)];
prob.bux =  inf * ones(2 * n + m + 3, 1);

prob.cones.type = [res.symbcon.MSK_CT_QUAD];
prob.cones.sub = [2, 3, zeros(1, m)];
for i = 1 : m
    prob.cones.sub(2 + i) = n + 3 + i;
end
prob.cones.subptr = [1];

[r,res]=mosekopt('minimize',prob);
x = res.sol.itr.xx(4 : n + 3, 1);
