function [ u,Balance_image ] = deblur_wav( image,A )
%%DEBLUR_WAV 此处显示有关此函数的摘要
%   此处显示详细说明
ATA = A'*A;
kappa = 0.5;
step_size_ista = 0.2;
step_size_ADMM = 1.618;
iter_time = 50;
lambda_ista = 0.8;


iter_time_ADMM = 50;
mu = 0.1 ;
lambda_ADMM = 0.01;

frame=1; % type of wavelet frame used: 0 is Haar; 1 is piecewise linear; 3 is piecewise cubic
Level=2; % level of decomposition, typically 1-4.
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel2D(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel2D(x,R,Level); % Frame reconstruction

len=length(D);

%%Analysis Approach
%% BY ADMM

u=image;
[m,n]=size(u);
d=W(image);
b=mul_wavcof(W(image),0,len);

for i = 1:iter_time_ADMM
    tmp = mu*WT(sub_wavcof(d,b,len));
    
    tmp = reshape(tmp,m*n,1);
    tmp = A'*tmp;
    tmp = reshape(tmp,m,n);
    
    u = (image+tmp);
    
    u = reshape(u,m*n,1);
    
    u = gmres(speye(m*n)*mu+ATA,u,10,1e-6,30);
    
    u = reshape(u,m,n);
    
    d = shrinkage(add_wavcof(W(u),b,len),lambda_ADMM,len);
    b = add_wavcof(b,mul_wavcof(sub_wavcof(W(u),d,len),step_size_ADMM,len),len);
end

%%Balanced Approach
%% By ISTA
g = W(image);

for i = 1:iter_time
    grad = sub_wavcof(add_wavcof(sub_wavcof(W(reshape(ATA*reshape(WT(g),m*n,1),m,n)),W(image),len),mul_wavcof(g,kappa,len),len),mul_wavcof(W(WT(g)),kappa,len),len);
    g = sub_wavcof(g,mul_wavcof(grad,step_size_ista,len),len);
    g = shrinkage(g,lambda_ista,len);
end

Balance_image = WT(g);

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

