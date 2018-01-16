addpath(genpath('util'))

%%
% Code for block gs method application in anisotropic poisson equation
% testing gridsize 64

for i=-4
    test(16,10^-i);
end

rmpath(genpath('util'))
