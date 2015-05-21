clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = true;
multiple_orbits = true;

%%function
[out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits);
                    
%%processing (plot/write to file)
figure('name','parameters over time')
t = out.tp;
subplot(5,2,1)
Rm = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
plot(t,Rm)
ylabel('$h$ $\left[m\right]$','interpreter','latex')
grid on
subplot(5,2,3)
Vm = sqrt(out.V(:,1).^2 + out.V(:,3).^2 + out.V(:,2).^2);
plot(t,Vm)
ylabel('$V$ $\left[\frac{m}{s}\right]$','interpreter','latex')
grid on
subplot(5,2,5)
hold on
am = sqrt((out.A(:,1) - out.Ag(:,1)).^2 + (out.A(:,2) - out.Ag(:,2)).^2 + (out.A(:,3) - out.Ag(:,3)).^2) / g_earth;
plot(t,am)
plot(xlim,[3,3],'-.','color','r');
plot(xlim,[control.a/g_earth,control.a/g_earth],'-.','color','g');
ylabel('$a_{human}$ $\left[\frac{m}{s^2}\right]$','interpreter','latex')
grid on
subplot(5,2,7)
plot(t,out.q)
ylabel('$\bar{q}$  $\left[Pa\right]$','interpreter','latex')
grid on
subplot(5,2,9)
hold on
plot(t,out.M)
plot(xlim,[5,5],'-.','color','r');
ylabel('$M$  $\left[-\right]$','interpreter','latex')
grid on
subplot(5,2,2)
plot(t,out.CL)
ylabel('$C_L$  $\left[-\right]$','interpreter','latex')
grid on
subplot(5,2,4)
plot(t,out.CD)
ylabel('$C_D$  $\left[-\right]$','interpreter','latex')
grid on
subplot(5,2,6)
plot(t,out.alpha*180/pi)
ylabel('$\alpha$  $\left[^circ\right]$','interpreter','latex')
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
axis([-(R_m + h_atm)*1.5 (R_m + h_atm)*1.5 -(R_m + h_atm)*1.5 (R_m + h_atm)*1.5])
plot(out.R(:,1),out.R(:,2))
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g')
theta_plot = out.theta0:0.00001:out.theta;
rk = out.a * (1- out.e^2) ./ (1 + out.e * cos(theta_plot));
polar(theta_plot+out.theta_p,rk,'k');
plot(out.rp*cos(out.theta_p),out.rp*sin(out.theta_p),'*');
plot(-out.ra*cos(out.theta_p),-out.ra*sin(out.theta_p),'d');
