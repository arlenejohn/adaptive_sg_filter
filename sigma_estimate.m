function [sigma] = sigma_estimate(t,type)
% the function estimates the signal standard deviation
% t is the median estimate of the signal
% type is the type of noise distribution ( type = 'G' for Gaussian, type = 'L' for laplacian, type = 'U' for
%   uniform )
if (type=='G')
    sigma = t/(0.6745*sqrt(2));
elseif (type=='L')
        sigma = (t/1.1461)*sqrt(2);
elseif (type=='U')
        sigma = t/((2-sqrt(2))*sqrt(3));
end