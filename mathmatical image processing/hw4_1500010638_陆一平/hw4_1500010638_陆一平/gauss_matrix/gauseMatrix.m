function [A] = gauseMatrix(m,n,gauss_rad,k)
%parameter


tmp = padarray(k, [m - (2 * gauss_rad + 1), n - (2 * gauss_rad + 1)], 'post');
tmp = circshift(tmp, [-gauss_rad, -gauss_rad]);
A = sparse(m * n, m * n);

for i = 1 : m
    A(i, :) = reshape(tmp, m * n, 1)';
    tmp = circshift(tmp, [1, 0]);
end

for i = 1 : gauss_rad
    A(i * m + 1 : (i + 1) * m, :) = circshift(A((i - 1) * m + 1 : i * m, :), [0, m]);
end



for i = gauss_rad + 1 : n - gauss_rad - 1
    A(i * m + 1: (i + 1) * m, (i - gauss_rad) * m + 1: (i + gauss_rad + 1) * m) = ...
        A(gauss_rad * m + 1: (gauss_rad + 1) * m, 1: (2 * gauss_rad + 1) * m);
end

for i = n - gauss_rad : n - 1
    A(i * m + 1 : (i + 1) * m, :) = circshift(A(1 : m, :), [0, m * i]);
end

% thanks to censhicong tells me how to calculate the matrix
end