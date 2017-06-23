function [ output_args ] = showseg( input_args )
%SHOWSEG show the segement of levelsetfunction

[m,n]=size(input_args);
output_args=zeros(m,n);
for i = 1:m
    for j = 1:n
        if input_args(i,j)>=0
            output_args(i,j)=256;
        end
    end
end



end

