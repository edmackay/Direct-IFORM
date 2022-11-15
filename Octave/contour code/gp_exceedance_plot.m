function gp_exceedance_plot(z,u,sigma,xi,conf)

% exceedance plot for GP model
%   z      - array of observations
%   u      - GP threshold parameter
%   sigma  - GP scale parameter
%   xi     - GP shape parameter

if nargin<5
    conf=0;
end

z=z(z>u);
n=length(z);

% empirical exceedance probabilities
z=sort(z);
k=1:n;
P=k/(n+1);

% GP exceedance probabilities
x=linspace(u,z(end),1000);
F=gp_cdf(x,xi,sigma,u);

% plot results
plot(z,1-P,'.','color',[1 1 1]*0.5)
hold on;
plot(x,1-F,'k')
set(gca,'yscale','log')
ylabel('Exceedance probability')
if conf==1
    % confidence bounds on exceedance probabilities
    a=k;
    b=n-k+1;
    Plow = betaincinv(0.025,a,b);
    Phigh = betaincinv(0.975,a,b);
    plot(z,1-Plow,'--','color',[1 1 1]*0.5)
    plot(z,1-Phigh,'--','color',[1 1 1]*0.5)
end
L=legend('Data','GP fit','location','southwest');
%L.Box='off';
set(L, 'Box','off');
