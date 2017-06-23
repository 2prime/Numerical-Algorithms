clear;
clc;
addpath('./profunc/');
% get test image
InImage = dir('./test_img/*.jpg');
num = size(InImage,1);

%denoise
for index = 1:num
    image=im2double(imread(['./test_img/' InImage(index).name]));
    sigma_2 = max(max(image)) / 100;
    image_noise = image +  7*sigma_2(1, 1) * randn(size(image,1),size(image,2),size(image,3));
    image_show=image;
    image_show(:,:,1)=denoise_wav(image_noise(:,:,1));
    image_show(:,:,2)=denoise_wav(image_noise(:,:,2));
    image_show(:,:,3)=denoise_wav(image_noise(:,:,3));
    
    image_show1=image;
    image_show1(:,:,1)=denoise_anl(image_noise(:,:,1));
    image_show1(:,:,2)=denoise_anl(image_noise(:,:,2));
    image_show1(:,:,3)=denoise_anl(image_noise(:,:,3));
    
    image_2 = denoise_TV(image_noise);
    saveas(imshow([image image_noise image_show image_show1 image_2]),['./result/wav_' InImage(index).name],'jpg')
end