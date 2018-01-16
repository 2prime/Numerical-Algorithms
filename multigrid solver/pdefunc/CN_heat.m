function [ output ] = CN_heat( init,n,t,tool )
%CN_HEAT 热方程CN格式
%   输入：init 初始值
%        n 时间网格数目
%        t  output_time

%%
time_step = t/n;
tmp_output = im_to_vec(init);

vec_lap = @(x) im_to_vec(lap_z(vec_to_im(x)));
nc_op = @(x) x-vec_lap(x)*time_step/2;

for i=1:n
    b = tmp_output + vec_lap(tmp_output)*time_step/2;
    if strcmp(tool.method,'gmres')
        tmp_output = gmres(nc_op,b,5);
    end
    
    if strcmp(tool.method,'pcg')
        tmp_output = pcg(nc_op,b,1e-8,100);
    end
    
    if strcmp(tool.method,'cg')
        tmp_output = cgs(nc_op,b,1e-8,100);
    end
    
    
end

output = vec_to_im(tmp_output);

end

