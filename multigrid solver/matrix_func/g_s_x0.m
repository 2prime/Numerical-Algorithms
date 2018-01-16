function [ x ] = g_s_x0( A,b,x0,end_tool )
%G_S_X0 gs iter start at x0 to solve Ax=b
%   end_tool: Iter_tool
%           iter_time: max iter time(the longest time to stop
%           iter_end: the ending condition

%%

D = diag(diag(A));
L = D - tril(A);
U = D - triu(A);
x = x0;

g = (D-L)\b;
g2 = (D-U)\b;


%%
for i = 1:end_tool.iter_time
    %iter step
    x =(D-L)\(U*x)+g;
    x =(D-U)\(L*x)+g2;

    if norm(A*x-b)/norm(b) < end_tool.iter_end
        fprintf('g-s method iter end at %g\n',end_tool.iter_end);
        break;
    end
end

fprintf('g-s error: %g\n',norm(A*x-b)/norm(b));


end

