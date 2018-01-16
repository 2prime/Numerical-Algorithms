function [ output ] = implicit_heat( init,n,t,tool)
%implicit_HEAT 热方程隐示格式
%   输入：init 初始值
%        n 时间网格数目
%        t  output_time

%%
time_step = t/n;
tmp_output = im_to_vec(init);

output = init;

grid_size = size(init,1)-1;

vec_lap = @(x) im_to_vec(lap_z(vec_to_im(x)));
im_op = @(x) x-vec_lap(x)*time_step;

tmp_z_output = im_to_vec(init(2:grid_size,2:grid_size));
lap_matrix = lapmatrix(grid_size).*time_step;
A = speye((grid_size-1)^2)-lap_matrix;

if strcmp(tool.method,'chol')
   R = full(chol(A)); 
end
for i=1:n
    if strcmp(tool.method,'chol')
        tmp_z_output = R'\(R\tmp_z_output);
    end
    
    if strcmp(tool.method,'test')
        tmp_z_output = pcg(A,tmp_z_output,1e-8,100);
    end
    
    if strcmp(tool.method,'gs')
        Iter_opt = Iter_tool;
        Iter_opt.iter_time = 1000;
        Iter_opt.iter_end = 1e-8;
        tmp_z_output = g_s(A,tmp_z_output,Iter_opt);
    end
    
    if strcmp(tool.method,'multigrid')
        output(2:grid_size,2:grid_size) = multigrid_heat(output(2:grid_size,2:grid_size),time_step,2);
    end
    
    if strcmp(tool.method,'gmres')
        tmp_output = gmres(im_op,tmp_output,5);
    end
    
    if strcmp(tool.method,'pcg')
        tmp_output = pcg(im_op,tmp_output,1e-8,100);
    end
    
    if strcmp(tool.method,'cg')
        tmp_output = cgs(im_op,tmp_output,1e-8,100);
    end
end

if strcmp(tool.method,'test')
    z_output = vec_to_im(tmp_z_output);
    output = zeros(grid_size+1);
    output(2:grid_size,2:grid_size) = z_output;
end

if strcmp(tool.method,'gs')
    z_output = vec_to_im(tmp_z_output);
    output = zeros(grid_size+1);
    output(2:grid_size,2:grid_size) = z_output;
end

if strcmp(tool.method,'chol')
    z_output = vec_to_im(tmp_z_output);
    output = zeros(grid_size+1);
    output(2:grid_size,2:grid_size) = z_output;
end

if strcmp(tool.method,'gmres')
    output = vec_to_im(tmp_output);
end

if strcmp(tool.method,'pcg')
    output = vec_to_im(tmp_output);
end

if strcmp(tool.method,'cg')
    output = vec_to_im(tmp_output);
end

end

