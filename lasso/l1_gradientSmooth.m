function [x, out] = l1_gradientSmooth(x0, A, b, mu, opts)

m = size(A, 1);
n = size(A, 2);

x = x0;

%calculate frequently used items
A_A = A' * A;
A_b = A' * b;

% some paramounts
mu0 = mu * 1e5;
smoothPara = 1e-4;
maxOutIter = 30;
maxInnerIter = 50;
precision = 1e-5;
muDescentRate = (1e5)^(-1/maxOutIter);
smoothDescentRate = 0.5;

% get err
ansNorm = norm(opts.ans);
err = zeros(1500, 1);
iterNum = 0;

for outIterCount = 1 : maxOutIter
    for innerIterCount = 1 : maxInnerIter
        
        % calculate \nabla g(x)
        g = A_A * x - A_b;
        for i = 1 : n
            g(i) = g(i) + mu0 * x(i)/sqrt(x(i)^2 + smoothPara);
        end
        
        step = 1e-3;
        val = func(x, A_A, A_b, mu0, smoothPara);
        while func(x - step * g, A_A, A_b, mu0, smoothPara) > val - step * 0.5 * g' * g
            step = step * 0.5;
        end
       
        % converge judgement
        if max(abs(g)) < precision
            break;
        end
                
        % calculate x - t \nabla g(x)
        x = x - step * g;
        
        iterNum = iterNum + 1;
        err(iterNum, 1) = norm(x - opts.ans) / (1 + ansNorm);
    end
    smoothPara = smoothPara * smoothDescentRate;
    mu0 = max(muDescentRate * mu0, mu);
end
out = struct('errMap', err);
end

function [value] = func(x, A_A, A_b, mu, smoothPara)
    value = 0.5 * x' * A_A * x - A_b' * x;
    for i = 1 : size(x, 1)
        value = value + mu * sqrt(x(i)^2 + smoothPara);
    end
end