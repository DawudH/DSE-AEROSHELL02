function [ g, p, T, rho, a ] = mars_atmosphere(h)
%[ g, p, T, rho, a ] = mars_atmosphere(h)
%   Generate atmosphere variables based on height and standard atmosphere
%   parameters as given in mars_atmosphere.m. Height can be a vector.
%   everything in SI units (m, s, kg, K)


% Load standard parameters
mars_standard_parameters

h_base = 1000*[4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,89,94,99,104];
T_base = [203.2,208.8,206.84,199.81,192.6,185.28,176.65,167.92,160.74,155.36,149.47,148.2,145.57,142.16,147.25,155.74,153.54,146.55,139.26,135.08,132.84];
p_base = [411.7,254.5,158.2,96.34,58.17,34.56,20.07,11.29,6.231,3.369,1.717,0.8915,0.4597,0.2307,0.1191,0.05869,0.0319,0.01736,0.009069,0.004631,0.002313];
rho_base = [0.01067,0.006416,0.004026,0.002538,0.00159,0.0009818,0.0005981,0.0003542,0.0002042,0.0001144,6.066e-05,3.178e-05,1.67e-05,8.607e-06,4.276e-06,2.104e-06,1.11e-06,6.244e-07,3.453e-07,1.815e-07,9.27e-08];


T = interp1(h_base, T_base, h, 'pchip', 'extrap');
p = interp1(h_base, p_base, h, 'pchip', 'extrap');
rho = interp1(h_base, rho_base, h, 'pchip', 'extrap');

%Heights outside the data region are given the value of the largest known
%height. Extrapolation yields dissatisfactory results.
T(h>h_base(end)) = T_base(end);
p(h>h_base(end)) = 0;
rho(h>h_base(end)) = 0;
T(h<h_base(1)) = T_base(1);
p(h<h_base(1)) = p_base(1);
rho(h<h_base(1)) = rho_base(1);

%Old version, not suitable for high altitudes
% T = T0 + dTdh.*h;
% p = p0*exp(p_grad.*h);
% rho = p./(Rm_mars*T);


g = G*M_mars./(h+r_mars).^2;
a = sqrt(gamma_mars*Rm_mars.*T);

end

