function [x, out] = l1_subGradient(x0, A, b, mu, opts)
x = x0;
m = size(A, 1);
n = size(A, 2);

%calculate frequently used items
A_A = A' * A;
A_b = A' * b;

threshold = 1e-3;
tmp1 = mu * ones(n, 1);
    
for iterCount = 1 : 2000   
    % calculate subgradient
    g = A_A * x - A_b;
    g(x > threshold) = g(x > threshold) + tmp1(x > threshold);
    g(x < threshold) = g(x < threshold) - tmp1(x < threshold);
    
    %backtrack search
%     step = 1;
%     rollback = 0.5;
%     wolfe = 0.01;
%     while func(x - step * g) > func(x) - wolfe * g' * (g * step)
%         step = step * rollback;
%     end
    
    % fixed step
%     step = 1e-2;   

    % fixed length
%     step = 1e-2 / norm(g, 2);

    % diminishing step size
    step = 1 / iterCount;
    
    %converge judgement
    x = x - g * step;
    if max(abs(g)) < 1e-5
        break
    end
end

