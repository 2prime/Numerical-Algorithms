clear;
clc;
addpath('./profunc/');
% get test image
InImage = dir('./testimg/*.jpg');
num = size(InImage,1);


%delur

InImage = dir('./test_img/*.jpg');
num=size(InImage,1);
for index = 2:2
    image=double(imread(['./test_img/' InImage(index).name]));
    [m,n,l]=size(image);
    [blured_image,result,anl,bal ]=blur(image);
    saveas(imshow(uint8([image blured_image anl bal result ])),['.\result\' InImage(index).name],'jpg')
    error = image-result;
    error_1 = reshape(error(:,:,1),m*n,1);error_2 = reshape(error(:,:,2),m*n,1);error_3 = reshape(error(:,:,3),m*n,1);
    image_1 = reshape(image(:,:,1),m*n,1);image_2 = reshape(image(:,:,1),m*n,1);image_3 = reshape(image(:,:,1),m*n,1);
    fprintf('Error rate:%f\n',(norm(error_1)+norm(error_2)+norm(error_3))/(norm(image_1)+norm(image_2)+norm(image_3)));
end