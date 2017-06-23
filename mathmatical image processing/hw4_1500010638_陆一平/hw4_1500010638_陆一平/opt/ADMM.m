function [ u ] = ADMM( f,W,lap,process,lambda,tol,max_iter_time )
fprintf('ADMM\n');
%ADMM for TV problem
learning_rate = 1.618;
mu = 0.1; %augmented
u = f;
d = W*u;
b=zeros(size(d));

size_tmp = size(W);
nabla_x = W(1:(size_tmp(1)/2),:);
nabla_y = W((size_tmp(1)/2)+1:(size_tmp(1)),:);

n=size(d,1);

%first step of ADMM
tmp = f+W'*(d-b).*mu;
k = size(W,2);
u=gmres(speye(k)+(W'*W).*mu,tmp,10,1e-6,30);
d=W*u+b;
gradx_tmp = nabla_x*u;
grady_tmp = nabla_y*u;
grad_tmp = gradx_tmp.*gradx_tmp + grady_tmp.*grady_tmp;
for i = 1:(n/2)
    if grad_tmp(i)>lambda/mu
        d(i)=d(i)*(1-(lambda/mu)/grad_tmp(i));
        d((n/2)+i)=d((n/2)+i)*(1-(lambda/mu)/grad_tmp(i));
    else
        d(i)=0;d((n/2)+i)=0;
    end
end
% d=sign(d).*max(abs(d)-lambda/mu,0);

% for i = 1:n
%    d(i)= Soft_threshold(d(i),lambda/mu);
% end

b=b+(W*u-d).*learning_rate;
tmp_matrix = process+(lap).*mu;

for i = 1:max_iter_time
    if norm(W*u-d)<tol*norm(f)
        fprintf('Compeled\n');
        break;
    end
    
    %homology
%     if mod(i,10) == 0
%         lambda=lambda*10;
%     end
    
    %ADMM iterate step
        tmp = f+W'*(d-b).*mu;
        u=gmres(tmp_matrix,tmp);
        d=W*u+b;
        gradx_tmp = nabla_x*u;
grady_tmp = nabla_y*u;
grad_tmp = gradx_tmp.*gradx_tmp + grady_tmp.*grady_tmp;
for i = 1:(n/2)
    if grad_tmp(i)>lambda/mu
        d(i)=d(i)*(1-(lambda/mu)/grad_tmp(i));
        d((n/2)+i)=d((n/2)+i)*(1-(lambda/mu)/grad_tmp(i));
    else
        d(i)=0;d((n/2)+i)=0;
    end
end
%         d=sign(d).*max(abs(d)-lambda/mu,0);
%         for i = 1:n
%             d(i)= Soft_threshold(d(i),lambda/mu);
%         end
        b=b+(W*u-d).*learning_rate;
    
end

end

function y = Soft_threshold(x,lambda)
if abs(x)<lambda
    y=0;
else
    if x>0
        y = x-lambda;
    else
        y = x+lambda;
    end
end
end

