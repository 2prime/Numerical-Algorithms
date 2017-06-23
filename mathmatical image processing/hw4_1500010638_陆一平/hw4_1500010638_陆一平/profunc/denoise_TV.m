function [image_show] = denoise(image)
    addpath('./gridfunc/');
    addpath('./opt/');


    [m,n]=size(image);

    Im1_vector = zeros(m*n,1);
    for i = 1:n
        Im1_vector((i-1)*m+1:i*m,:)= image(:,i);
    end
    
    fprintf('begin processing\n');

    % solve total variation model by ADMM
    W=[gradfuncx(m,n); gradfuncy(m,n)];
    lap = laplacfunc(m,n);

    u = ADMM(Im1_vector,W,lap,speye(m*n),0.001,0.0005,20);

    % image for showing
    image_show = zeros(size(image));
    for i = 1:n
        image_show(:,i) = u((i-1)*m+1:i*m,:);
    end
    

end
  



