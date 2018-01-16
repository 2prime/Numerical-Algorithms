function [ tmp ] = lap_z( output )
%LAP_Z: laplace with zero padding
%   laplace = u_xx+u_yy
%   five point finite difference scheme
%%
lap = laplace(output);
tmp = zeros(size(lap,1)+2);
tmp(2:size(lap,1)+1,2:size(lap,1)+1)=lap;
end

