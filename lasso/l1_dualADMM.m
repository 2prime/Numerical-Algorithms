function [x, out] = l1_dualADMM(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);


%calculate frequently used items
A_A = A' * A;
AA_ = A * A';
A_b = A' * b;

% some paramounts
mu0 = 1e2;
maxIter = 50;
precision = 1e-5;
muDescentRate = (1e2 / mu)^(-1/(maxIter - 1));

if (mu > mu0)
    mu0 = mu;
    muDescentRate = 1;
end


% useful
g2 = mu * ones(n, 1);

% get errMap
ansNorm = norm(opts.ans);
err = zeros(1000, 1);
iterNum = 0;

beta = 1e-2;
z = zeros(n,1);
y = -x0 / beta;
%y_old=y;
%p=3;
for iterCount = 1 : maxIter
    %y=(p/(iterCount+p))*y_old+(iterCount/(iterCount+p))*y;
    [L, D] = ldl(eye(m) + AA_ * beta);
    x = L'\(D\(L\(-b + beta * A *(y - z))));
    tmp = y - A' * x;
    z = sign(tmp) .* min(abs(tmp), mu0);
    %y_old=y;
    y = y - 1.618 * (A' * x + z);
    
    % calculate err
    iterNum = iterNum + 1;
    err(iterNum, 1) = norm(y * beta - opts.ans)/(1e-6 + ansNorm);
    
    mu0 = max(mu0 * muDescentRate, mu);
end
x = y * beta;
out = struct('errMap', err);