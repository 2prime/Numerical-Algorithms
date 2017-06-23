function [ u ] = denoise_anl( image )
%DENOISE_ANL 此处显示有关此函数的摘要
%   此处显示详细说明
iter_time = 50;
mu = 0.1 ;
lambda = 0.1;
step_size = 1.618;


frame=1; % type of wavelet frame used: 0 is Haar; 1 is piecewise linear; 3 is piecewise cubic
Level=2; % level of decomposition, typically 1-4.
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel2D(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel2D(x,R,Level); % Frame reconstruction

len=length(D);
u=image;
d=W(image);
b=mul_wavcof(W(image),0,len);

for i = 1:iter_time
    u = (image+mu*WT(sub_wavcof(d,b,len)))/(1+mu);
    d = shrinkage(add_wavcof(W(u),b,len),lambda,len);
    b = add_wavcof(b,mul_wavcof(sub_wavcof(W(u),d,len),step_size,len),len);
end


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


