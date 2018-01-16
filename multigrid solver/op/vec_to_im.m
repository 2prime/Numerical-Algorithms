function [ im ] = vec_to_im( vec )
%VEC_TO_IM the inverse of im_to_vec

n = sqrt(size(vec,1));
im = zeros(n,n);

for i = 1:n
    im(:,i) = vec(1+(i-1)*n:i*n);
end

end

