function [ vec ] = im_to_vec( im )
%IM_TO_VEC 

n = size(im,1);
vec = zeros(n^2,1);

for i = 1:n
    vec(1+(i-1)*n:i*n)=im(:,i);
end

end

