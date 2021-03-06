clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = false;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

%%function
[out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
                    
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
ylabel('$V$ $\left[m s^{-1}\right]$','interpreter','latex')
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
grid on
hold on
plot(t,out.theta,'color',cc(2,:))
ylabel('$\theta$  $\left[deg\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,5)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.gamma,'color',cc(2,:))
ylabel('$\gamma$  $\left[deg\right]$','interpreter','latex')
%xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end

subplot(3,3,9)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.alpha*180/pi,'color',cc(2,:))
ylabel('$\alpha$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end


subplot(3,3,8)
xlim([min(t) max(t)])
grid on
hold on
plot(t,out.phi*180/pi,'color',cc(2,:))
ylabel('$\mu$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
if isfield(out,'tkep') 
    plot([out.tkep, out.tkep],ylim,'-.','color',cc(1,:),'LineWidth',1.4); 
end
if (export_figures)
    matlab2tikz('.\LaTeX\orbit_results2.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end


% plot orbit
% circle plot:

figure('name','Orbit')
%grid on
axis equal
hold on
%axis([-(R_m + h_atm)*2 (R_m + h_atm)*2 -(R_m + h_atm)*2 (R_m + h_atm)*2])
axis([-(R_m + h_atm)*1.2 (R_m + h_atm)*1.2 -(R_m + h_atm)*1.2 (R_m + h_atm)*1.2])
if hypkep
    % plot the hyperbolic part
    theta_plot = out.theta0:0.001:out.theta;
    rk = out.a * (1- out.e^2) ./ (1 + out.e * cos(theta_plot));
    h1 = polar(theta_plot+out.theta_p,rk); 
    set(h1,'color',cc(4,:))
    legend_str{1} = 'Hyperbolic kepler';
end
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
        if hypkep
            legend([h1, L1, h2, L2, h3, h4],legend_str,'location','east');
        else
            legend([L1, h2, L2, h3, h4],legend_str(2:end),'location','east');
        end
    else
        if hypkep
            legend([h1, L1, h3, h4],legend_str);
        else
            legend([L1, h3, h4],legend_str(2:end));
        end
    end
    % plot location of start landing phase:
    point1 = 1.0e+06 * [-3.274898658418613, -1.907769574406631, 0];
    point = 1.0e+06 * [-3.393522184766882,  -0.278472149680601,  0];
    plot(point(1),point(2),'x','color','k','markers',12)
    plot(point1(1),point1(2),'x','color','k','markers',12)
    if (export_figures)
        matlab2tikz('.\LaTeX\orbit2.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
    end
