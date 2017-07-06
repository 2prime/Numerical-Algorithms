function [x, out] = l1_linearADMM(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);

% some paramounts
mu0 = mu * 1e5;
maxIter = 100;
muDescentRate = (1e5)^(-1/(maxIter - 5));


% get errMap
ansNorm = norm(opts.ans);
err = zeros(1000, 1);
iterNum = 0;

beta = 3e-2;
x = x0;
z = A * x;
y = beta * ones(m, 1);

for iterCount = 1 : maxIter
    gradx = beta * A' * (A * x - z - y);
    for i = 1 : n
        gradx(i) = gradx(i) + mu0 * x(i)/sqrt(x(i)^2 + 1 / (iterCount^2));
    end
    
    step = 1;
    fx = func(x, A, y, z, mu0, beta);
    while func(x - step * gradx, A, y, z, mu0, beta) > fx - step * gradx' * gradx * 0.5
        step = step * 0.5;
    end
    
    x = x - step * gradx;
    z = (b + beta * (A * x - y)) / sqrt(1 + beta);
    y = y - 1.618 * (A * x - z);
    
    
    % calculate err
    iterNum = iterNum + 1;
    err(iterNum, 1) = norm(x - opts.ans)/(1e-6 + ansNorm);
    
    mu0 = max(mu0 * muDescentRate, mu);
end
out = struct('errMap', err);
end

function [value] = func(x, A, y, z, mu, beta)
    value = mu * sum(abs(x)) + 0.5 * beta * sum((A * x - z - y).^2);
end