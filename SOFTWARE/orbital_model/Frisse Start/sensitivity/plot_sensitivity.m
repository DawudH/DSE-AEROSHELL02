close all
clear
clc

addpath('../')
addpath('../../../matlab2tikz')
constants 

export_figures = false;

 out_1   = load('aerocapture_rho_1.mat');
 out_1_1 = load('aerocapture_rho_1_1.mat');
 out_0_9 = load('aerocapture_rho_0_9.mat');


%% processing for the aerocapture
figure('name','parameters over time for the aerocapture')
t_1 = out_1.out.tp(1:end-1);
t_1_1 = out_1_1.out.tp(1:end-1);
t_0_9 = out_0_9.out.tp(1:end-1);

cc = parula(7);
c1 = cc(1,:);
c2 = cc(3,:);
c3 = cc(5,:);

xlimits = [0 max(t_1_1)];

if export_figures
    index1 = 1:50:length(t_1);
    index2 = 1:50:length(t_1_1);
    index3 = 1:50:length(t_0_9);
else
    index1 = 1:1:length(t_1);
    index2 = 1:1:length(t_1_1);
    index3 = 1:1:length(t_0_9);
end


markerspace1 = 3:500:length(t_1)-3;
markerspace2 = 3:500:length(t_1_1)-3;
markerspace3 = 3:500:length(t_0_9)-3;

linespec = {'*','o','s'};

subplot(3,3,1)
xlim(xlimits)
grid on
hold on
Rm_1 = sqrt(out_1.out.R((1:end-1),1).^2 + out_1.out.R((1:end-1),3).^2 + out_1.out.R((1:end-1),2).^2) - R_m;
plot(t_1(index1),Rm_1(index1),'color',c1)
plot(t_1(markerspace1),Rm_1(markerspace1),linespec{1},'color',c1)
Rm_1_1 = sqrt(out_1_1.out.R((1:end-1),1).^2 + out_1_1.out.R((1:end-1),3).^2 + out_1_1.out.R((1:end-1),2).^2) - R_m;
plot(t_1_1(index2),Rm_1_1(index2),'color',c2)
plot(t_1_1(markerspace2),Rm_1_1(markerspace2),linespec{2},'color',c2)
Rm_0_9 = sqrt(out_0_9.out.R((1:end-1),1).^2 + out_0_9.out.R((1:end-1),3).^2 + out_0_9.out.R((1:end-1),2).^2) - R_m;
plot(t_0_9(index3),Rm_0_9(index3),'color',c3)
plot(t_0_9(markerspace3),Rm_0_9(markerspace3),linespec{3},'color',c3)
ylabel('$h$ $\left[m\right]$','interpreter','latex')

subplot(3,3,4)
xlim(xlimits)
grid on
hold on
Vm = sqrt(out_1.out.V(:,1).^2 + out_1.out.V(:,3).^2 + out_1.out.V(:,2).^2);
plot(t_1(index1),Vm((index1)),'color',c1)
plot(t_1(markerspace1),Vm(markerspace1),linespec{1},'color',c1)
Vm = sqrt(out_1_1.out.V(:,1).^2 + out_1_1.out.V(:,3).^2 + out_1_1.out.V(:,2).^2);
plot(t_1_1(index2),Vm(index2),'color',c2)
plot(t_1_1(markerspace2),Vm(markerspace2),linespec{2},'color',c2)
Vm = sqrt(out_0_9.out.V(:,1).^2 + out_0_9.out.V(:,3).^2 + out_0_9.out.V(:,2).^2);
plot(t_0_9(index3),Vm(index3),'color',c3)
plot(t_0_9(markerspace3),Vm(markerspace3),linespec{3},'color',c3)

subplot(3,3,7)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.a_human_mag((index1))/g_earth,'color',c1)
plot(t_1_1(index2),out_1_1.out.a_human_mag(index2)/g_earth,'color',c2)
plot(t_0_9(index3),out_0_9.out.a_human_mag(index3)/g_earth,'color',c3)
plot(t_1(markerspace1),out_1.out.a_human_mag(markerspace1)/g_earth,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.a_human_mag(markerspace2)/g_earth,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.a_human_mag(markerspace3)/g_earth,linespec{3},'color',c3)
ylabel('$a_{astronaut}$ $\left[g_e\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
plot(xlim,[3,3],'-.','color',cc(5,:),'LineWidth',1.4);

subplot(3,3,3)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.q((index1)),'color',c1)
plot(t_1_1(index2),out_1_1.out.q(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.q(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.q(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.q(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.q(markerspace3),linespec{3},'color',c3)
ylabel('$\bar{q}$  $\left[Pa\right]$','interpreter','latex')


subplot(3,3,6)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.M((index1)),'color',c1)
plot(t_1_1(index2),out_1_1.out.M(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.M(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.M(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.M(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.M(markerspace3),linespec{3},'color',c3)
ylabel('$M$  $\left[-\right]$','interpreter','latex')
plot(xlim,[5,5],'-.','color',cc(5,:),'LineWidth',1.4);


subplot(3,3,2)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.theta(index1),'color',c1)
plot(t_1_1(index2),out_1_1.out.theta(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.theta(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.theta(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.theta(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.theta(markerspace3),linespec{3},'color',c3)
ylabel('$\theta$  $\left[deg\right]$','interpreter','latex')
h = legend('Nominal trajectory','Trajectory with 10% more density','Trajectory with 10% less density','location','north','orientation','horizontal');
h.Position = [0.243115519562129 0.955832322216614 0.533872688020616 0.0198986978848316];

subplot(3,3,5)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.gamma(index1),'color',c1)
plot(t_1_1(index2),out_1_1.out.gamma(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.gamma(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.gamma(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.gamma(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.gamma(markerspace3),linespec{3},'color',c3)
ylabel('$\gamma$  $\left[deg\right]$','interpreter','latex')

subplot(3,3,9)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.alpha(index1)*180/pi,'color',c1)
plot(t_1_1(index2),out_1_1.out.alpha(index2)*180/pi,'color',c2)
plot(t_0_9(index3),out_0_9.out.alpha(index3)*180/pi,'color',c3)
plot(t_1(markerspace1),out_1.out.alpha(markerspace1)*180/pi,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.alpha(markerspace2)*180/pi,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.alpha(markerspace3)*180/pi,linespec{3},'color',c3)
ylabel('$\alpha$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')


subplot(3,3,8)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.phi(index1)*180/pi,'color',c1)
plot(t_1_1(index2),out_1_1.out.phi(index2)*180/pi,'color',c2)
plot(t_0_9(index3),out_0_9.out.phi(index3)*180/pi,'color',c3)
plot(t_1(markerspace1),out_1.out.phi(markerspace1)*180/pi,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.phi(markerspace2)*180/pi,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.phi(markerspace3)*180/pi,linespec{3},'color',c3)
ylabel('$\mu$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')

if (export_figures)
    matlab2tikz('.\orbit_sensitivity_capture.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end



%% processing (plot/write to file)

out_1_c = out_1;
out_1_1_c = out_1_1;
out_0_9_c = out_0_9;

out_1   = load('entry_rho_1.mat');
out_1_1 = load('entry_rho_1_1.mat');
out_0_9 = load('entry_rho_0_9.mat');


figure('name','parameters over time for the entry')
t_1 = out_1.out.tp(1:end-1);
t_1_1 = out_1_1.out.tp(1:end-1);
t_0_9 = out_0_9.out.tp(1:end-1);

xlimits = [0 max(t_1_1)];

if export_figures
    index1 = 1:50:length(t_1);
    index2 = 1:50:length(t_1_1);
    index3 = 1:50:length(t_0_9);
else
    index1 = 1:1:length(t_1);
    index2 = 1:1:length(t_1_1);
    index3 = 1:1:length(t_0_9);
end


markerspace1 = 3:500:length(t_1)-3;
markerspace2 = 3:500:length(t_1_1)-3;
markerspace3 = 3:500:length(t_0_9)-3;

linespec = {'*','o','s'};

subplot(3,3,1)
xlim(xlimits)
grid on
hold on
Rm_1 = sqrt(out_1.out.R((1:end-1),1).^2 + out_1.out.R((1:end-1),3).^2 + out_1.out.R((1:end-1),2).^2) - R_m;
plot(t_1(index1),Rm_1(index1),'color',c1)
plot(t_1(markerspace1),Rm_1(markerspace1),linespec{1},'color',c1)
Rm_1_1 = sqrt(out_1_1.out.R((1:end-1),1).^2 + out_1_1.out.R((1:end-1),3).^2 + out_1_1.out.R((1:end-1),2).^2) - R_m;
plot(t_1_1(index2),Rm_1_1(index2),'color',c2)
plot(t_1_1(markerspace2),Rm_1_1(markerspace2),linespec{2},'color',c2)
Rm_0_9 = sqrt(out_0_9.out.R((1:end-1),1).^2 + out_0_9.out.R((1:end-1),3).^2 + out_0_9.out.R((1:end-1),2).^2) - R_m;
plot(t_0_9(index3),Rm_0_9(index3),'color',c3)
plot(t_0_9(markerspace3),Rm_0_9(markerspace3),linespec{3},'color',c3)
ylabel('$h$ $\left[m\right]$','interpreter','latex')

subplot(3,3,4)
xlim(xlimits)
grid on
hold on
Vm = sqrt(out_1.out.V(:,1).^2 + out_1.out.V(:,3).^2 + out_1.out.V(:,2).^2);
plot(t_1(index1),Vm((index1)),'color',c1)
plot(t_1(markerspace1),Vm(markerspace1),linespec{1},'color',c1)
Vm = sqrt(out_1_1.out.V(:,1).^2 + out_1_1.out.V(:,3).^2 + out_1_1.out.V(:,2).^2);
plot(t_1_1(index2),Vm(index2),'color',c2)
plot(t_1_1(markerspace2),Vm(markerspace2),linespec{2},'color',c2)
Vm = sqrt(out_0_9.out.V(:,1).^2 + out_0_9.out.V(:,3).^2 + out_0_9.out.V(:,2).^2);
plot(t_0_9(index3),Vm(index3),'color',c3)
plot(t_0_9(markerspace3),Vm(markerspace3),linespec{3},'color',c3)

subplot(3,3,7)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.a_human_mag((index1))/g_earth,'color',c1)
plot(t_1_1(index2),out_1_1.out.a_human_mag(index2)/g_earth,'color',c2)
plot(t_0_9(index3),out_0_9.out.a_human_mag(index3)/g_earth,'color',c3)
plot(t_1(markerspace1),out_1.out.a_human_mag(markerspace1)/g_earth,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.a_human_mag(markerspace2)/g_earth,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.a_human_mag(markerspace3)/g_earth,linespec{3},'color',c3)
ylabel('$a_{astronaut}$ $\left[g_e\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')
plot(xlim,[3,3],'-.','color',cc(5,:),'LineWidth',1.4);

subplot(3,3,3)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.q((index1)),'color',c1)
plot(t_1_1(index2),out_1_1.out.q(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.q(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.q(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.q(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.q(markerspace3),linespec{3},'color',c3)
ylabel('$\bar{q}$  $\left[Pa\right]$','interpreter','latex')


subplot(3,3,6)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.M((index1)),'color',c1)
plot(t_1_1(index2),out_1_1.out.M(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.M(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.M(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.M(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.M(markerspace3),linespec{3},'color',c3)
ylabel('$M$  $\left[-\right]$','interpreter','latex')
plot(xlim,[5,5],'-.','color',cc(5,:),'LineWidth',1.4);


subplot(3,3,2)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.theta(index1),'color',c1)
plot(t_1_1(index2),out_1_1.out.theta(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.theta(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.theta(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.theta(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.theta(markerspace3),linespec{3},'color',c3)
ylabel('$\theta$  $\left[deg\right]$','interpreter','latex')
h = legend('Nominal trajectory','Trajectory with 10% more density','Trajectory with 10% less density','location','north','orientation','horizontal');
h.Position = [0.243115519562129 0.955832322216614 0.533872688020616 0.0198986978848316];

subplot(3,3,5)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.gamma(index1),'color',c1)
plot(t_1_1(index2),out_1_1.out.gamma(index2),'color',c2)
plot(t_0_9(index3),out_0_9.out.gamma(index3),'color',c3)
plot(t_1(markerspace1),out_1.out.gamma(markerspace1),linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.gamma(markerspace2),linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.gamma(markerspace3),linespec{3},'color',c3)
ylabel('$\gamma$  $\left[deg\right]$','interpreter','latex')

subplot(3,3,9)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.alpha(index1)*180/pi,'color',c1)
plot(t_1_1(index2),out_1_1.out.alpha(index2)*180/pi,'color',c2)
plot(t_0_9(index3),out_0_9.out.alpha(index3)*180/pi,'color',c3)
plot(t_1(markerspace1),out_1.out.alpha(markerspace1)*180/pi,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.alpha(markerspace2)*180/pi,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.alpha(markerspace3)*180/pi,linespec{3},'color',c3)
ylabel('$\alpha$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')


subplot(3,3,8)
xlim(xlimits)
grid on
hold on
plot(t_1(index1),out_1.out.phi(index1)*180/pi,'color',c1)
plot(t_1_1(index2),out_1_1.out.phi(index2)*180/pi,'color',c2)
plot(t_0_9(index3),out_0_9.out.phi(index3)*180/pi,'color',c3)
plot(t_1(markerspace1),out_1.out.phi(markerspace1)*180/pi,linespec{1},'color',c1)
plot(t_1_1(markerspace2),out_1_1.out.phi(markerspace2)*180/pi,linespec{2},'color',c2)
plot(t_0_9(markerspace3),out_0_9.out.phi(markerspace3)*180/pi,linespec{3},'color',c3)
ylabel('$\mu$  $\left[deg\right]$','interpreter','latex')
xlabel('$t$ $\left[s\right]$','interpreter','latex')

if (export_figures)
    matlab2tikz('.\orbit_sensitivity_entry.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end



%% Plot entry orbit at mars
figure('name','Orbit')
axis equal
hold on
axis([-(R_m + h_atm)*1.2 (R_m + h_atm)*1.2 -(R_m + h_atm)*1.2 (R_m + h_atm)*1.2])

% planet and atmosphere
    theta_plot = 0:0.01:2*pi;
    radius_mars = ones(1,length(theta_plot)) * R_m;
    radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
    h3 = polar(theta_plot,radius_mars); 
    set(h3,'color',cc(6,:))
    %legend_str{end+1} = 'Mars mean-equatorial radius';
    h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
    set(h4,'color',cc(5,:))
    %legend_str{end+1} = 'Mars atmosphere limit';
    set(gca,'Visible','off')
    set(gcf,'color',[1 1 1])
    
    

% if isfield(out,'tkep')
%     ccc = autumn(5);
%     %first atmos section
%     index = find(t == out.tkep);
%     L1 = plot(out.R(1:10:index,1),out.R(1:10:index,2),'color',cc(1,:));
%     legend_str{2} = 'First pass through atmosphere';
%     plot(out.R(1:150:index,1),out.R(1:150:index,2),'v','color',cc(1,:));
%     % plot kepler eliptical part
%     theta_plot = out.okep.theta:0.02:(2*pi-out.okep.theta);
%     rk = out.okep.a * (1- out.okep.e^2) ./ (1 + out.okep.e .* cos(theta_plot));
%     h2 = polar(theta_plot+out.okep.theta_p,rk); 
%     set(h2,'color',cc(3,:))
%     legend_str{3} = 'Elliptical kepler';
%     % plot second atmos section
%     L2 = plot(out.R(index+1:10:end,1),out.R(index+1:10:end,2),'color',ccc(2,:)); 
%     plot(out.R(index+1:150:end,1),out.R(index+1:150:end,2),'d','color',ccc(2,:)); 
%     legend_str{4} = 'Second pass through atmosphere';
% else
%     % just plot the atmos part
%     L1 = plot(out.R(:,1),out.R(:,2),'color',cc(1,:),'LineWidth',1.2);
%     legend_str{2} = 'Pass through atmosphere';
% end
% 
%     if exist('L2','var')
%         if hypkep
%             legend([h1, L1, h2, L2, h3, h4],legend_str,'location','east');
%         else
%             legend([L1, h2, L2, h3, h4],legend_str(2:end),'location','east');
%         end
%     else
%         if hypkep
%             legend([h1, L1, h3, h4],legend_str);
%         else
%             legend([L1, h3, h4],legend_str(2:end));
%         end
%     end
%     % plot location of start landing phase:
%     point1 = 1.0e+06 * [-3.274898658418613, -1.907769574406631, 0];
%     point = 1.0e+06 * [-3.393522184766882,  -0.278472149680601,  0];;
%     plot(point(1),point(2),'x','color','k','markers',12)
%     plot(point1(1),point1(2),'x','color','k','markers',12)
%     if (export_figures)
%         matlab2tikz('.\orbit_sensitivity_entry_mars.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
%     end