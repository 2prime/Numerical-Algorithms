function [x, out] = l1_gurobi(x0, A, b, mu, opts)
clear model

x = x0;

m = size(A, 1);
n = size(A, 2);

model.Q = sparse(diag([zeros(1, n), 0.5 * ones(1, m), zeros(1, n)]));
model.obj = [zeros(1, m + n), mu * ones(1, n)];
model.modelsense = 'min';

model.A = zeros(m + 2 * n, m + 2 * n);
model.A(1 : m, 1 : m + n) = [A, -1 * eye(m)];
for i = 1 : n
    model.A(m + i, i) = 1;
    model.A(m + i, m + n + i) = 1;
    model.A(m + n + i, i) = -1;
    model.A(m + n + i, m + n + i) = 1;
end
model.A = sparse(model.A);
model.rhs = [b; zeros(2 * n, 1)];
model.sense = char([61 * ones(1, m), 62 * ones(1, 2 * n)]);

result = gurobi(model);
x = result.x(1 : n, 1);