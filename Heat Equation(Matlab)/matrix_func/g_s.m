function [ x ] = g_s( A,b,end_tool )
%G_S Gauss-Seidel Iter To solve Ax=b
%   end_tool: Iter_tool
%           iter_time: max iter time(the longest time to stop
%           iter_end: the ending condition

%%

D = diag(diag(A));
L = D - tril(A);
U = D - triu(A);
x = zeros(size(A,1),1);

g = (D-L)\b;

%%
for i = 1:end_tool.iter_time
    %iter step
    x =(D-L)\(U*x)+g;

    if norm(A*x-b)/norm(b) < end_tool.iter_end
        fprintf('g-s method iter end at %g\n',end_tool.iter_end);
        break;
    end
end

fprintf('g-s error: %g\n',norm(A*x-b)/norm(b));


end

