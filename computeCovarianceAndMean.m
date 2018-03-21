function [c,mu] = computeCovarianceAndMean(y,sigma_w,r)
%% Compute c and mu 
c  = (1/sigma_w^2+1./r(:)).^(-1);                                            % diagonal-elements of cov. matrix 
inv= ifft2(y);                                                               % fourier-inverse of the observation
mu = (1/sigma_w^2)*c.*inv(:);                                                % mean of the posterior-distribution
