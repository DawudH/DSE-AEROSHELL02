clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = true;

%%function
[out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits);
                    
%% processing (plot/write to file)
figure('name','parameters over time')
t = out.tp;

cc = parula(7);

subplot(3,3,1)
xlim([min(t) max(t)])
grid on
hold on
Rm = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
plot(t,Rm,'color',cc(2,:))
ylabel('$h$ $\left[m\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4);
end


subplot(3,3,4)
xlim([min(t) max(t)])
grid on
hold on
Vm = sqrt(out.V(:,1).^2 + out.V(:,3).^2 + out.V(:,2).^2);
plot(t,Vm,'color',cc(2,:))
ylabel('$V$ $\left[\frac{m}{s}\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end


subplot(3,3,7)
xlim([min(t) max(t)])
grid on
hold on
am = sqrt((out.A(:,1) - out.Ag(:,1)).^2 + (out.A(:,2) - out.Ag(:,2)).^2 + (out.A(:,3) - out.Ag(:,3)).^2) / g_earth;
plot(t,am,'color',cc(2,:))
ylabel('$a_{astronaut}$ $\left[g_e\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
plot(xlim,[3,3],'-.','color',cc(5,:),'LineWidth',1.4);
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,3)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.q,'color',cc(2,:))
ylabel('$\bar{q}$  $\left[Pa\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,6)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.M,'color',cc(2,:))
ylabel('$M$  $\left[-\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
plot(xlim,[5,5],'-.','color',cc(5,:),'LineWidth',1.4);
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,2)
xlim([min(t) max(t)])
ylim([min(out.CL)-0.005*abs(min(out.CL)) max(out.CL)+0.005*abs(max(out.CL))])
grid on
hold on
plot(t,out.CL,'color',cc(2,:))
ylabel('$C_L$  $\left[-\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,5)
xlim([min(t) max(t)])
ylim([min(out.CD)-0.001*abs(min(out.CD)) max(out.CD)+0.001*abs(max(out.CD))])
grid on
hold on
plot(t,out.CD,'color',cc(2,:))
ylabel('$C_D$  $\left[-\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,8)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.alpha*180/pi,'color',cc(2,:))
ylabel('$\alpha$  $\left[^\circ\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end


subplot(3,3,9)
xlim([min(t) max(t)])
ylim([-abs(control.dalphadt)*180/pi abs(control.dalphadt)*180/pi])
ylim([-0.05 0.32])
grid on
hold on
plot(t,out.dAlpha_dt*180/pi,'color',cc(2,:))
ylabel('$\frac{d \alpha}{dt}$  $\left[\frac{^\circ}{s}\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end
matlab2tikz('.\LaTeX\orbit_results2.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);


% plot orbit
% circle plot:

figure('name','Orbit')
%grid on
axis equal
hold on
%axis([-(R_m + h_atm)*2 (R_m + h_atm)*2 -(R_m + h_atm)*2 (R_m + h_atm)*2])
axis([-(R_m + h_atm)*1.2 (R_m + h_atm)*1.2 -(R_m + h_atm)*1.1 (R_m + h_atm)*1.1])

% plot the hyperbolic part
theta_plot = out.theta0:0.001:out.theta;
rk = out.a * (1- out.e^2) ./ (1 + out.e * cos(theta_plot));
h1 = polar(theta_plot+out.theta_p,rk); 
set(h1,'color',cc(4,:))
legend_str{1} = 'Hyperbolic kepler';
if isfield(out,'tkep')
    ccc = autumn(5);
    %first atmos section
    index = find(t == out.tkep);
    L1 = plot(out.R(1:10:index,1),out.R(1:10:index,2),'color',cc(1,:));
    legend_str{2} = 'First pass through atmosphere';
    plot(out.R(1:150:index,1),out.R(1:150:index,2),'v','color',cc(1,:));
    % plot kepler eliptical part
    theta_plot = out.okep.theta:0.02:(2*pi-out.okep.theta);
    rk = out.okep.a * (1- out.okep.e^2) ./ (1 + out.okep.e .* cos(theta_plot));
    h2 = polar(theta_plot+out.okep.theta_p,rk); 
    set(h2,'color',cc(3,:))
    legend_str{3} = 'Elliptical kepler';
    % plot second atmos section
    L2 = plot(out.R(index+1:10:end,1),out.R(index+1:10:end,2),'color',ccc(2,:)); 
    plot(out.R(index+1:150:end,1),out.R(index+1:150:end,2),'d','color',ccc(2,:)); 
    legend_str{4} = 'Second pass through atmosphere';
else
    % just plot the atmos part
    L1 = plot(out.R(:,1),out.R(:,2),'color',cc(1,:),'LineWidth',1.2);
    legend_str{2} = 'Pass through atmosphere';
end
% planet and atmosphere
    theta_plot = 0:0.01:2*pi;
    radius_mars = ones(1,length(theta_plot)) * R_m;
    radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
    h3 = polar(theta_plot,radius_mars); 
    set(h3,'color',cc(6,:))
    legend_str{end+1} = 'Mars mean-equatorial radius';
    h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
    set(h4,'color',cc(5,:))
    legend_str{end+1} = 'Mars atmosphere limit';
    set(gca,'Visible','off')
    set(gcf,'color',[1 1 1])
    if exist('L2','var')
        legend([h1, L1, h2, L2, h3, h4],legend_str,'location','east');
    else
        legend([h1, L1, h3, h4],legend_str);
    end
    matlab2tikz('.\LaTeX\orbit2.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);