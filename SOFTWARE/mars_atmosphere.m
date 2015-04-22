function [ g, p, T, rho, a ] = mars_atmosphere(h)
%[ g, p, T, rho, a ] = mars_atmosphere(h)
%   Generate atmosphere variables based on height and standard atmosphere
%   parameters as given in mars_atmosphere.m. Height can be a vector.


% Load standard parameters
mars_standard_parameters

T = T0 + dTdh.*h;
a = sqrt(gamma_mars*Rm_mars.*T);
p = p0*exp(p_grad.*h);
rho = p./(Rm_mars*T);
g = GM_mars./(h+r_mars).^2;

end

