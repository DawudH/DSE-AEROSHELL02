clc
clear all
close all

%%Input constants & variables
variables
h = 200e3;
dv2 = -8;

% booleans
use_control = false;
multiple_orbits = true;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

%%function
%[out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
load('../../Aero model/orbits/orbit_iteration_1_1.mat')
mu = G*M_mars;
ra = out.okep.a*(1-out.okep.e^2)/(1-out.okep.e);
v = sqrt(mu*(2/ra-1/out.okep.a));
rp = R_m + h;
a_n = (rp+ra)/2;
e_n = 1-rp/a_n;
T_n = 2*pi * sqrt(a_n^3 / mu)/3600
v_n = sqrt(mu*(2/ra-1/a_n));
dv1 = v_n - v
r_e = R_m + h_atm;
v_e = v_n + dv2;
a_e = (2/ra-v_e^2/mu)^-1
e_e = ra/a_e-1
v_entry = sqrt(mu*(2/r_e-1/a_e))
theta_entry = acos((a_e*(1-e_e^2)-r_e)/(r_e*e_e))
theta0 = theta_entry + out.okep.theta_p
b_e = a_e*sqrt(1-e_e^2);
x = r_e*cos(theta0);
rc_e = b_e^2*(a_e*e_e+x)/(a_e^2*sqrt(b_e^2-b_e^2*(a_e*e_e+x)^2/a_e^2))

%gamma = 17.2/180*pi;
%theta = (2*pi-out.okep.theta);
% theta_p = acos((cos(theta0)*cos(gamma)-sin(theta0)*sin(gamma)+rc_e*(sin(gamma)*cos(theta0)+cos(gamma)*sin(theta0)))/sqrt(1+rc_e^2))
% theta = theta0-theta_p
%e_e = -(1/2)*(cos(theta)*r0+sqrt(cos(theta)^2*r0^2-2*cos(theta)*r0^2+4*cos(theta)*r0*ra+r0^2+4*r0*ra-4*ra^2)+r0)/(cos(theta)*r0-ra)
%e_e = (1/2)*(-cos(theta)*r0-r0+sqrt(cos(theta)^2*r0^2-2*cos(theta)*r0^2+4*cos(theta)*r0*ra+r0^2+4*r0*ra-4*ra^2))/(cos(theta)*r0-ra)
%v_entry = sqrt(mu*(2/r_e-1/a_e));
%a_e = (2/r_e-v_entry^2/mu)^-1
%e_e = ra/a_e-1
%v_e = sqrt(mu*(2/ra-1/a_e));
%b_e = a_e*sqrt(1-e_e^2);
%V_unit = [1,rc_e,0]/norm([1,rc_e,0]);
%V_unit = (rotz(-theta_p)*V_unit')'
%theta_e = acos( (a_e*(1-e_e^2)-r) / (e_e*r) );
ccc = jet(8);
figure('name','eliptic orbit')
axis equal
hold on
%axis([-(R_m + h_atm)*2 (R_m + h_atm)*2 -(R_m + h_atm)*2 (R_m + h_atm)*2])
%axis([-(R_m + h_atm)*1.2 (R_m + h_atm)*1.2 -(R_m + h_atm)*1.2 (R_m + h_atm)*1.2])
theta_plot = out.okep.theta:0.02:(pi);
rk = out.okep.a * (1- out.okep.e^2) ./ (1 + out.okep.e .* cos(theta_plot));
h2 = polar(theta_plot+out.okep.theta_p,rk); 
set(h2,'color',ccc(3,:))
%legend_str{3} = 'Elliptical kepler';
theta_plot = 0:0.02:2*pi;
rk = a_n * (1- e_n^2) ./ (1 + e_n .* cos(theta_plot));
hn = polar(theta_plot+out.okep.theta_p,rk); 
set(hn,'color',ccc(1,:))
theta_plot = pi:0.02:2*pi;
rk = a_e * (1- e_e^2) ./ (1 + e_e .* cos(theta_plot));
he = polar(theta_plot+out.okep.theta_p,rk); 
set(he,'color',ccc(4,:))
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
h3 = polar(theta_plot,radius_mars); 
%set(h3,'color',ccc(6,:))
%legend_str{end+1} = 'Mars mean-equatorial radius';
h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
%set(h4,'color',cc(5,:))
%legend_str{end+1} = 'Mars atmosphere limit';
set(gca,'Visible','off')
set(gcf,'color',[1 1 1])