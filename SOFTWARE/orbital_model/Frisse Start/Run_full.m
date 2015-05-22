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
                    
%% processing (plot/write to file)
figure('name','parameters over time')
t = out.tp;

cc = parula(7);

subplot(5,2,1)
xlim([min(t) max(t)])
grid on
hold on
Rm = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
plot(t,Rm,'color',cc(2,:))
ylabel('$h$ $\left[m\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4);
end


subplot(5,2,3)
xlim([min(t) max(t)])
grid on
hold on
Vm = sqrt(out.V(:,1).^2 + out.V(:,3).^2 + out.V(:,2).^2);
plot(t,Vm,'color',cc(2,:))
ylabel('$V$ $\left[\frac{m}{s}\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end


subplot(5,2,5)
xlim([min(t) max(t)])
grid on
hold on
am = sqrt((out.A(:,1) - out.Ag(:,1)).^2 + (out.A(:,2) - out.Ag(:,2)).^2 + (out.A(:,3) - out.Ag(:,3)).^2) / g_earth;
plot(t,am,'color',cc(2,:))
ylabel('$a_{human}$ $\left[\frac{m}{s^2}\right]$','interpreter','latex')
plot(xlim,[3,3],'-.','color',cc(5,:),'LineWidth',1.4);
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(5,2,7)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.q,'color',cc(2,:))
ylabel('$\bar{q}$  $\left[Pa\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(5,2,9)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.M,'color',cc(2,:))
ylabel('$M$  $\left[-\right]$','interpreter','latex')
plot(xlim,[5,5],'-.','color',cc(5,:),'LineWidth',1.4);
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(5,2,2)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.CL,'color',cc(2,:))
ylabel('$C_L$  $\left[-\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(5,2,4)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.CD,'color',cc(2,:))
ylabel('$C_D$  $\left[-\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(5,2,6)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.alpha*180/pi,'color',cc(2,:))
ylabel('$\alpha$  $\left[^circ\right]$','interpreter','latex')
if exist('out.tkep','var') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end



% plot orbit
% circle plot:
figure('name','Orbit')
grid on
axis equal
hold on
axis([-(R_m + h_atm)*2 (R_m + h_atm)*2 -(R_m + h_atm)*2 (R_m + h_atm)*2])
% plot the hyperbolic part
theta_plot = out.theta0:0.00001:out.theta;
rk = out.a * (1- out.e^2) ./ (1 + out.e * cos(theta_plot));
h1 = polar(theta_plot+out.theta_p,rk); 
set(h1,'color',cc(4,:))

theta_plot = 0:0.01:2*pi;
if isfield(out,'tkep')
    ccc = autumn(5);
    %first atmos section
    index = find(t == out.tkep);
    plot(out.R(1:index,1),out.R(1:index,2),'color',ccc(2,:),'LineWidth',1.1); 
    % plot kepler eliptical part
    rk = out.okep.a * (1- out.okep.e^2) ./ (1 + out.okep.e .* cos(theta_plot));
    h2 = polar(theta_plot+out.okep.theta_p,rk); 
    set(h2,'color',cc(4,:))
    % plot second atmos section
    plot(out.R(index+1:end,1),out.R(index+1:end,2),'color',ccc(2,:),'LineWidth',1.1); 
else
    % just plot the atmos part
    plot(out.R(:,1),out.R(:,2))
end
% planet and atmosphere
    theta_plot = 0:0.01:2*pi;
    radius_mars = ones(1,length(theta_plot)) * R_m;
    radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
    h3 = polar(theta_plot,radius_mars); 
    set(h3,'color',cc(6,:))
    h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
    set(h4,'color',cc(5,:))
