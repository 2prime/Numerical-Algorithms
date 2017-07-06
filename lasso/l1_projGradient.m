function [x, out] = l1_projGradient(x0, A, b, mu, opts)
m = size(A, 1);
n = size(A, 2);

%calculate frequently used items
A_A = A' * A;
A_b = A' * b;

x1 = zeros(n, 1);
x2 = zeros(n, 1);
x1(x0 > 0) = x0(x0 > 0);
x2(x0 < 0) = -x0(x0 < 0);

x0 = zeros(n, 1);

for iterCount = 1 : 10000
   % calculate grad
   g1 = A_A * (x1 - x2) - A_b + mu * ones(n, 1);
   g2 = A_A * (x2 - x1) + A_b + mu * ones(n, 1);
   
   %converge judgement
   if max(abs(g1)) < 1e-5 && max(abs(g2)) < 1e-5
       break;
   end
   
   % calculate step
   step = (sum(g1.^2) + sum(g2.^2)) / (g1' * A_A * g1 + g2' * A_A * g2 - 2 * g1' * A_A * g2 + 1e-10);
   
   % get next point
   x1 = x1 - step * g1;
   x2 = x2 - step * g2;
   
   % projection
   x1(x1 < 0) = x0(x1 < 0);
   x2(x2 < 0) = x0(x2 < 0);
end
x = x1 - x2;