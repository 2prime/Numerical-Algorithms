function [x, out] = l1_proxGradient(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);

x = x0;

%calculate frequently used items
A_A = A' * A;
A_b = A' * b;

% fixed / first step size for BB
step = 1 / norm(A_A);

% some paramounts
mu0 = mu * 1e5;
threshold = 1e-5;
maxOutIter = 16;
maxInnerIter = 50;
precision = 1e-5;
muDescentRate = (1e5)^(-1/maxOutIter);

% useful
g2 = mu * ones(n, 1);

% get errMap
ansNorm = norm(opts.ans);
err = zeros(1500, 1);
iterNum = 0;

for outIterCount = 1 : maxOutIter
    softThreshold = step * mu0;
    shrinkVect = softThreshold * ones(n, 1);
    for innerIterCount = 1 : maxInnerIter
        % calculate \nabla g(x)
        g = A_A * x - A_b;
        
        % BB
%         if iterCount > 1
%             s = x - prev_x;
%             step = (s' * s) / (s' * (g - prev_g));
%         end
%         prev_g = g;
%         prev_x = x;
        
        % calculate x - t \nabla g(x)
        x_tmp = x - step * g;
        
        % converge judgement
        g(x > threshold) = g(x > threshold) + g2(x > threshold);
        g(x < threshold) = g(x < threshold) - g2(x < threshold);
        if max(abs(g)) < precision
            break;
        end
        
        % execute prox operator
        x = zeros(n, 1);
        x(x_tmp > softThreshold) = x_tmp(x_tmp > softThreshold) - shrinkVect(x_tmp > softThreshold);
        x(x_tmp < -softThreshold) = x_tmp(x_tmp < -softThreshold) + shrinkVect(x_tmp < -softThreshold);
        
        % calculate err
        iterNum = iterNum + 1;
        err(iterNum, 1) = norm(x - opts.ans)/(1 + ansNorm);
    end
    if (mu0 == mu)
        break;
    end
    mu0 = max(muDescentRate * mu0, mu);
end
out = struct('errMap', err);