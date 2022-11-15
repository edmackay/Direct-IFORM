function F = gp_cdf(x,xi,sigma,threshold)
% Generalized Pareto cumulative density function (cdf).
% xi = shape parameter (xi>0 -> heavy tail, xi<0 -> bounded tail)
% sigma = scale parameter

F = zeros(size(x));
z = x-threshold;

% calculate CDF for interval where x>threshold
ind = z > 0;
if xi ~=0.0
    F(ind) = 1-(1+xi*z(ind)/sigma).^(-1/xi);
else
    F(ind) = 1-exp(-z(ind)/sigma);
end

% set CDF to zero for interval where x<threshold
F(~ind) = 0;

% set CDF to one where x exceeds upper bound
if xi<0
    ind=z>-sigma/xi;
    F(ind)=1;
end