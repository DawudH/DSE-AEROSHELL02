clc
clear
close all

%%Input constants & variables
variables
%escape velocity at boundary of atmosphere
V_esc = sqrt(2*G*M_mars/(R_m+h_atm));

%%function
[out] = full_orbit(R, V, V_esc, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth);
                    
%%processing (plot/write to file)
figure('name','parameters over time')
t = out.tp;
subplot(5,1,1)
Rm = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2);
plot(t,Rm)
grid on
subplot(5,1,2)
Vm = sqrt(out.V(:,1).^2 + out.V(:,3).^2 + out.V(:,2).^2);
plot(t,Vm)
grid on
subplot(5,1,3)
am = sqrt((out.A(:,1) - out.Ag(:,1)).^2 + (out.A(:,2) - out.Ag(:,2)).^2 + (out.A(:,3) - out.Ag(:,3)).^2) / g_earth;
plot(t,am)
grid on
subplot(5,1,4)
plot(t,out.q)
grid on
subplot(5,1,5)
hold on
plot(t,out.M)
plot(xlim,[5,5],'-.','color','r');
grid on


% plot orbit
% circle plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
figure('name','Orbit')
grid on
axis equal
hold on
plot(out.R(:,1),out.R(:,2))
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g')