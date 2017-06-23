clear;
addpath './img'
addpath './levelsetfunc'
addpath './denoise'
InImage = dir('./img/*.jpg');
num = size(InImage,1);

%% levelset method
for index =1:3
    image=double(imread(['./img/' InImage(index).name]));
    image=image(:,:,1);
    [v,g]=levelset(image,20000,20);
    
    %% show image
    v_ = showseg(v);
    v=showlev(v);
    [m,n]=size(v);
    image_new = zeros(m,n,3);
    v_new = zeros(m,n,3);
    image_new(:,:,1)=image.*(ones(m,n)-v);
    image_new(:,:,2)=image.*(ones(m,n)-v)+256*v;
    image_new(:,:,3)=image.*(ones(m,n)-v);
    v_new(:,:,1)=v_;v_new(:,:,2)=v_;v_new(:,:,3)=v_;
    
    saveas(imshow([uint8(image_new) uint8(v_new)]) ,['./test/20000' InImage(index).name],'jpg');

end