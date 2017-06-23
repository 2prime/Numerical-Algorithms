function [ levelset ] = showlev( image )
%SHOWLEV show the levelset of the image
%authorrized by @2pimre
%2017.3


[m,n]=size(image);
levelset=zeros(m,n);
for i = 2:m-1
    for j = 2:n-1
        if image(i,j)*image(i-1,j)<=0||image(i,j)*image(i+1,j)<=0||image(i,j)*image(i,j+1)<=0||image(i,j)*image(i,j-1)<=0
            levelset(i,j)=1;
        end
        %if abs(image(i,j))<5
        %    levelset(i,j)=1;
        %end
    end
end

end

