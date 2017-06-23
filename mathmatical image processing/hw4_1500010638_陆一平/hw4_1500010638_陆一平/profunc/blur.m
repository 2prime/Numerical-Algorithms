function [blured_image,result,anl,balance]=blur(image)
gauss_rad = 7;
sigma_1 = 2;
k = fspecial('Gaussian', [2 * gauss_rad + 1, 2 * gauss_rad + 1], sigma_1);
[m,n,l]=size(image);
addpath('./gridfunc/');
addpath('./opt/');
addpath('./gauss_matrix/');

sigma_2 = max(max(image)) / 100;

for i =  1 : l
    tmp = padarray(image(:, :, i), [gauss_rad, gauss_rad], 'circular', 'both');
    blured_image(:, :, i) = conv2(tmp, k, 'valid') + sigma_2(1, 1, i) * randn(m, n, 1);
end


W=[gradfuncx(m,n); gradfuncy(m,n)];
lap = laplacfunc(m,n);
process_matrix=gauseMatrix(m,n,gauss_rad,k);

anl = image;
balance = image;

for my_index = 1:l
    [anl(:,:,my_index),balance(:,:,my_index)]=deblur_wav(blured_image(:,:,my_index),process_matrix);
    
    Im1_vector = reshape(blured_image(:, :, my_index), m*n, 1);
    
    fprintf('begin processing\n');

    % solve total variation model by ADMM
    
    u = ADMM(Im1_vector,W,lap,process_matrix,3,0.001,100);

    % image for showing
    result(:, :, my_index) = reshape(u, m, n);
end
    

end
