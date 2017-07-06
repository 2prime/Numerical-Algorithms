function [x, out] = l1_augmentedLagrangian(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);

x = x0;

% some paramounts
mu0 = 1e2;
maxOutIter = 12;
maxInnerIter = 4;
muDescentRate = (1e2 / mu)^(-1/maxOutIter);

if (mu > mu0)
    mu0 = mu;
    muDescentRate = 1;
end

% get errMap
ansNorm = norm(opts.ans);
err = zeros(1000, 1);
iterNum = 0;

z = zeros(m, 1);
eta = 1e-2;

for outIterCount = 1 : maxOutIter
    for innerIterCount = 1 : maxInnerIter
        % Newton method
        for i = 1 : 3
        
            % Calculate grad
            tmp = eta * A' * z - x;
            grad = z + b + A * softThreshold(tmp, eta * mu0);

            % Calculate Hessian
            A0 = A(:, abs(tmp) > eta * mu0);
            hessian = eye(m) + eta * A0 * A0';

            % LDL solves Newton equation
%             [L,D] = ldl(hessian);
%             y = L'\(D\(L\grad));
            y = hessian \ grad;
        
            % Backtracking
            step = 1;
            hz = h(z, A, b, x, eta, mu0);
            while h(z - step * y, A, b, x, eta, mu0) > hz - step * y' * grad * 0.1
                step = step * 0.5;
            end
            z = z - step * y;
        end
        x = softThreshold(x - eta * A' * z, eta * mu0);
        
        % calculate err
        iterNum = iterNum + 1;
        err(iterNum, 1) = norm(x - opts.ans)/(1e-6 + ansNorm);
    end
    if (mu0 == mu)
        break;
    end
    mu0 = max(muDescentRate * mu0, mu);
end
out = struct('errMap', err);
end

function [x] = softThreshold(x0, mu)
    x = sign(x0) .* max(abs(x0) - mu, 0);
end

function [value] = h(z, A, b, x, eta, mu)
    tmp = softThreshold(A' * z - x / eta, mu);
    value = 0.5 * z' * z + b' * z + eta * 0.5 * tmp' * tmp;
end