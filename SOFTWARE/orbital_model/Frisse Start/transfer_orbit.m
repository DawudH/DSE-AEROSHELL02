clear
close all
clc

au = 1.496e11; %m
r_e = 1.013 * au; %m
r_m = 1.53 * au; %m

G = 6.67384*10^-11;
M_sun = 1.988435e30; % kg

V_e = sqrt(G * M_sun / r_e);
V_m = sqrt(G * M_sun / r_m);

delta_V_hohmann = sqrt( G*M_sun * (2/r_e - 2 / (r_m + r_e))) - V_e;
delta_V = delta_V_hohmann;
V_sc_p = delta_V + V_e;

a =  (2 / r_e - V_sc_p^2 / ( G * M_sun))^(-1);
e = 1- r_e / a;
b = sqrt(a^2 - (e*a)^2);

A_tot = pi * a * b

V_sc_m = sqrt( G * M_sun * (2/r_m - 1/a));
disp(['Velocity at mars: ' num2str(V_sc_m - V_m)])

theta = acos( (a * (1-e^2) / r_m - 1) / e);

fun = @(TH) 1/2 * (a * (1-e^2) ./ (1 + e*cos(TH))).^2;
A = integral(fun,0,theta)

%A = 2*sqrt(e^2 - 1) * a * atanh( (e - 1) * tan(theta/2) / sqrt(e^2-1))
t = 2*A / sqrt(a*G*M_sun*(1-e^2)) / 3600 / 24

figure('name','transfer orbit')
hold on
axis equal
theta_plot = 0:0.01:2*pi;
r_sc = a * (1-e^2) ./ (1 + e * cos(theta_plot) );
r_mars = ones(1,length(theta_plot)) * r_m;
r_earth = ones(1,length(theta_plot)) * r_e;
polar(theta_plot,r_mars); 
polar(theta_plot,r_earth); 
polar(theta_plot,r_sc); 