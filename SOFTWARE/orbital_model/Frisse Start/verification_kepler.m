clc
%%initialise
clear
close all

variables
atm.rhomesh = zeros(size(atm.rhomesh));

%%Hyperbolic kepler
[out] = hyperbolic_kepler(R,V,A,G,M_mars,R_m,h_atm,dt_kep_init);

%% numerical simulation
a_prev = A;
i = 1;
CL = 0;
CD = 0;
while R(i,2) > -(R_m+h_atm)
    [out_a] = in_atmosphere( V(i,:), R(i,:), A(i,:), a_prev, J(i,:), atm, CL, CD, dt_atmos, R_m, Omega_m, S, m );
    a_prev = A(i,:);
    R(i+1,:) = out_a.R;
    V(i+1,:) = out_a.V;
    A(i+1,:) = out_a.A;
    J(i+1,:) = out_a.J;
    i = i+1;
end
% circle plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
%radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
figure('name','Orbit')
grid on
axis equal
hold on
axis([-(R_m + h_atm)*1.5 (R_m + h_atm)*1.5 -(R_m + h_atm)*1.5 (R_m + h_atm)*1.5])
plot(R(:,1),R(:,2))
polar(theta_plot,radius_mars,'r');
%polar(theta_plot,radius_mars_atmos,'g')
theta_plot = -out.param.theta0:0.00001:-out.end.theta+0.1022*pi;
rk = out.param.a * (1- out.param.e^2) ./ (1 + out.param.e * cos(theta_plot));
polar(theta_plot+out.param.theta_p,rk,'k');
