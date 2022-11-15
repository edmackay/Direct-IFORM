function [sigma, xi]=gp_fit(z,min_xi,max_xi)

N=length(z);
zmax=max(z);

% first guess is moment estimator of sigma and xi
m=mean(z);
v=var(z);
xi0=0.5*(1-m^2/v);
if xi0<min_xi
    xi0=min_xi;
elseif xi0>max_xi
    xi0=max_xi;
end
sigma0=m*(1-xi0);
if xi0<0 && zmax>-sigma0/xi0
    sigma0=-xi0*zmax;
end
p0=[xi0, sigma0];

% maximise likelihood starting from first guess
opts=optimset('display','off');
params=fminsearch(@(p)GP_negloglike(p,z),p0,opts);
xi=params(1);
sigma=params(2);

% likelihood function
    function NLOGL = GP_negloglike(PARAMS,z)
        Xi = PARAMS(1);
        Sig = PARAMS(2);
        
        if Sig<0 || Xi<min_xi || Xi>max_xi || (Xi<0 && zmax>-Sig/Xi) 
            NLOGL = Inf;
        elseif abs(Xi) > eps
            NLOGL = N*log(Sig) + (1+1/Xi).*sum(log1p(Xi*z/Sig));
        else
            NLOGL = N*log(Sig) + sum(z/Sig);
        end
    end
end
