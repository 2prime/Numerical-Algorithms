function [ output_image ] = denoise_wav( image )
%DENOISE: use to denoise a picture 

addpath('./2D');
kappa = 0.5;
step_size = 0.01;
iter_time = 40;
lambda = 0.0008;

frame=1; % type of wavelet frame used: 0 is Haar; 1 is piecewise linear; 3 is piecewise cubic
Level=2; % level of decomposition, typically 1-4.
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel2D(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel2D(x,R,Level); % Frame reconstruction

len=length(D);
g = W(image);

for i = 1:iter_time
    grad = sub_wavcof(add_wavcof(sub_wavcof(W(WT(g)),W(image),len),mul_wavcof(g,kappa,len),len),mul_wavcof(W(WT(g)),kappa,len),len);
    g = sub_wavcof(g,mul_wavcof(grad,step_size,len),len);
    g = shrinkage(g,lambda,len);
end

output_image = WT(g);


end

function [wavcof3] = add_wavcof(wavcof1,wavcof2,len)
wavcof3 = wavcof2;
for kn=1:len-1
    for km=1:len-1
        if kn>1 || km>1
            wavcof3{1}{kn,km}=wavcof1{1}{kn,km}+wavcof2{1}{kn,km};
        end
    end
end
for kn=1:len-1
    for km=1:len-1
        if kn==1 && km==1
            wavcof3{2}{kn,km}=wavcof1{2}{kn,km}+wavcof2{2}{kn,km};
        end
        if kn>1 || km>1
            wavcof3{2}{kn,km}=wavcof1{2}{kn,km}+wavcof2{2}{kn,km};
        end
    end
end
end


function [wavcof3] = sub_wavcof(wavcof1,wavcof2,len)
wavcof3 = wavcof2;
for kn=1:len-1
    for km=1:len-1
        if kn>1 || km>1
            wavcof3{1}{kn,km}=wavcof1{1}{kn,km}-wavcof2{1}{kn,km};
        end
    end
end
for kn=1:len-1
    for km=1:len-1
        if kn==1 && km==1
            wavcof3{2}{kn,km}=wavcof1{2}{kn,km}-wavcof2{2}{kn,km};
        end
        if kn>1 || km>1
            wavcof3{2}{kn,km}=wavcof1{2}{kn,km}-wavcof2{2}{kn,km};
        end
    end
end
end


function [wavcof3] = mul_wavcof(wavcof1,k,len)
wavcof3 = wavcof1;
for kn=1:len-1
    for km=1:len-1
        if kn>1 || km>1
            wavcof3{1}{kn,km}=k*wavcof1{1}{kn,km};
        end
    end
end

for kn=1:len-1
    for km=1:len-1
        if kn==1 && km==1
            wavcof3{2}{kn,km}=k*wavcof1{2}{kn,km};
        end
        if kn>1 || km>1
            wavcof3{2}{kn,km}=k*wavcof1{2}{kn,km};
        end
    end
end
end

function [wavcof3] = shrinkage(wavcof1,lambda,len)
wavcof3 = wavcof1;
for kn=1:len-1
    for km=1:len-1
        if kn>1 || km>1
            wavcof3{1}{kn,km}= max(abs(wavcof1{1}{kn,km})-lambda,0).*sign(wavcof1{1}{kn,km});
        end
    end
end

for kn=1:len-1
    for km=1:len-1
        if kn==1 && km==1
            wavcof3{2}{kn,km}= max(abs(wavcof1{2}{kn,km})-lambda,0).*sign(wavcof1{2}{kn,km});
        end
        if kn>1 || km>1
            wavcof3{2}{kn,km}= max(abs(wavcof1{2}{kn,km})-lambda,0).*sign(wavcof1{2}{kn,km});
        end
    end
end
end
