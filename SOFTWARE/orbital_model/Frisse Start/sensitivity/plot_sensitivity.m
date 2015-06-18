close all
clear
clc

addpath('../')
addpath('../../../matlab2tikz')
constants 

export_figures = true;

 out_1   = load('aerocapture_rho_1.mat');
 out_1_1 = load('aerocapture_rho_1_1.mat');
 out_0_9 = load('aerocapture_rho_0_9.mat');

c_m = [161 37 27]/255;

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
Vm = sqrt(out_1.out.V(:,1).^2 + out_1.out.V(:,3).^2 + out_1.out.V(:,2).^2) /1000;
plot(t_1(index1),Vm((index1)),'color',c1)
plot(t_1(markerspace1),Vm(markerspace1),linespec{1},'color',c1)
Vm = sqrt(out_1_1.out.V(:,1).^2 + out_1_1.out.V(:,3).^2 + out_1_1.out.V(:,2).^2) /1000;
plot(t_1_1(index2),Vm(index2),'color',c2)
plot(t_1_1(markerspace2),Vm(markerspace2),linespec{2},'color',c2)
Vm = sqrt(out_0_9.out.V(:,1).^2 + out_0_9.out.V(:,3).^2 + out_0_9.out.V(:,2).^2) /1000;
plot(t_0_9(index3),Vm(index3),'color',c3)
plot(t_0_9(markerspace3),Vm(markerspace3),linespec{3},'color',c3)
ylabel('$V$ $\left[km s^{-1}\right]$','interpreter','latex')

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
ylabel('$a_{aero}$ $\left[g_e\right]$','interpreter','latex')
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
t_1 = out_1.out.tp;
t_1_1 = out_1_1.out.tp;
t_0_9 = out_0_9.out.tp;

xlimits = [0 max(t_1_1)];

if export_figures
    index1 = 1:50:length(t_1);
    index2 = 1:50:length(t_1_1);
    index3 = 1:50:length(t_0_9);
    index1o = 1:10:length(t_1);
    index2o = 1:10:length(t_1_1);
    index3o = 1:10:length(t_0_9);
else
    index1 = 1:1:length(t_1);
    index2 = 1:1:length(t_1_1);
    index3 = 1:1:length(t_0_9);
    index1o = 1:1:length(t_1);
    index2o = 1:1:length(t_1_1);
    index3o = 1:1:length(t_0_9);
end


markerspace1 = 3:500:length(t_1)-3;
markerspace2 = 3:500:length(t_1_1)-3;
markerspace3 = 3:500:length(t_0_9)-3;

linespec = {'*','o','s'};

subplot(3,3,1)
xlim(xlimits)
grid on
hold on
Rm_1 = sqrt(out_1.out.R(:,1).^2 + out_1.out.R(:,3).^2 + out_1.out.R(:,2).^2) - R_m;
plot(t_1(index1),Rm_1(index1),'color',c1)
plot(t_1(markerspace1),Rm_1(markerspace1),linespec{1},'color',c1)
Rm_1_1 = sqrt(out_1_1.out.R(:,1).^2 + out_1_1.out.R(:,3).^2 + out_1_1.out.R(:,2).^2) - R_m;
plot(t_1_1(index2),Rm_1_1(index2),'color',c2)
plot(t_1_1(markerspace2),Rm_1_1(markerspace2),linespec{2},'color',c2)
Rm_0_9 = sqrt(out_0_9.out.R(:,1).^2 + out_0_9.out.R(:,3).^2 + out_0_9.out.R(:,2).^2) - R_m;
plot(t_0_9(index3),Rm_0_9(index3),'color',c3)
plot(t_0_9(markerspace3),Rm_0_9(markerspace3),linespec{3},'color',c3)
ylabel('$h$ $\left[m\right]$','interpreter','latex')

subplot(3,3,4)
xlim(xlimits)
grid on
hold on
Vm = sqrt(out_1.out.V(:,1).^2 + out_1.out.V(:,3).^2 + out_1.out.V(:,2).^2) /1000;
plot(t_1(index1),Vm((index1)),'color',c1)
plot(t_1(markerspace1),Vm(markerspace1),linespec{1},'color',c1)
Vm = sqrt(out_1_1.out.V(:,1).^2 + out_1_1.out.V(:,3).^2 + out_1_1.out.V(:,2).^2) /1000;
plot(t_1_1(index2),Vm(index2),'color',c2)
plot(t_1_1(markerspace2),Vm(markerspace2),linespec{2},'color',c2)
Vm = sqrt(out_0_9.out.V(:,1).^2 + out_0_9.out.V(:,3).^2 + out_0_9.out.V(:,2).^2) /1000;
plot(t_0_9(index3),Vm(index3),'color',c3)
plot(t_0_9(markerspace3),Vm(markerspace3),linespec{3},'color',c3)
ylabel('$V$ $\left[km s^{-1}\right]$','interpreter','latex')

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
ylabel('$a_{aero}$ $\left[g_e\right]$','interpreter','latex')
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
figure('name','Entry orbit')
axis equal
hold on
%

% planet and atmosphere
    theta_plot = 0:0.01:2*pi;
    radius_mars = ones(1,length(theta_plot)) * R_m;
    radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
    h3 = polar(theta_plot,radius_mars); 
    set(h3,'color',c_m)
    %legend_str{end+1} = 'Mars mean-equatorial radius';
    h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
    set(h4,'color',cc(5,:))
    %legend_str{end+1} = 'Mars atmosphere limit';
    set(gca,'Visible','off')
    set(gcf,'color',[1 1 1])
 
% plot trajcetories
L1 = plot(out_1.out.R(index1o,1),out_1.out.R(index1o,2),'color',c1);
L2 = plot(out_1_1.out.R(index2o,1),out_1_1.out.R(index2o,2),'color',c2);
L3 = plot(out_0_9.out.R(index3o,1),out_0_9.out.R(index3o,2),'color',c3);
    % markers
    L1m = plot(out_1.out.R(markerspace1,1),out_1.out.R(markerspace1,2),'*','color',c1);
    L2m = plot(out_1_1.out.R(markerspace2,1),out_1_1.out.R(markerspace2,2),'o','color',c2);
    L3m = plot(out_0_9.out.R(markerspace3,1),out_0_9.out.R(markerspace3,2),'d','color',c3);

% determine the three landing locations
point1 = out_1.out.R(end,:);
dt2 = (out_1_1_c.out.tp(end) + out_1_1.out.tp(end)) - (out_1_c.out.tp(end) + out_1.out.tp(end));
dt3 = (out_0_9_c.out.tp(end) + out_0_9.out.tp(end)) - (out_1_c.out.tp(end) + out_1.out.tp(end));
point2 = rotz(dt2*omega_m*180/pi)*point1';
point3 = rotz(dt3*omega_m*180/pi)*point1';

L4 = plot(point1(1),point1(2),'x','color',c1,'markers',12);
L5 = plot(point2(1),point2(2),'x','color',c2,'markers',12);
L6 = plot(point3(1),point3(2),'x','color',c3,'markers',12);
rotate(L1,[0, 0, 1],-75);
rotate(L2,[0, 0, 1],-75);
rotate(L3,[0, 0, 1],-75);
rotate(L1m,[0, 0, 1],-75);
rotate(L2m,[0, 0, 1],-75);
rotate(L3m,[0, 0, 1],-75);
rotate(L4,[0, 0, 1],-75);
rotate(L5,[0, 0, 1],-75);
rotate(L6,[0, 0, 1],-75);
axis([-(R_m + h_atm)*0.7 (R_m + h_atm)*0.7 (R_m + h_atm)*0.75 (R_m + h_atm)*1])
legend([L1, L2, L3, h3, h4],'Nominal trajectory','Trajectory with 10% more density','Trajectory with 10% less density','Surface of Mars', 'Boundary of the atmosphere','location','south','orientation','horizontal')
if (export_figures)
    matlab2tikz('.\orbit_sensitivity_entry_mars.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end



%% Plot the complete orbit! :D
out_p   = load('parking.mat');
t1 = out_1_c.out.tp(1:end-1); % aerocapture
t2 = out_1.out.tp; % entry

if export_figures
    index1 = 1:50:length(t1);
    index2 = 1:50:length(t2);
else
    index1 = 1:1:length(t1);
    index2 = 1:1:length(t2);
end
figure('name','Complete orbit')
axis equal
hold on
%axis([-(R_m + h_atm)*1.2 (R_m + h_atm)*1.2 -(R_m + h_atm)*1.2 (R_m + h_atm)*1.2])

% planet and atmosphere
    theta_plot = 0:0.01:2*pi;
    radius_mars = ones(1,length(theta_plot)) * R_m;
    radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
    h3 = polar(theta_plot,radius_mars); 
    set(h3,'color',c_m)
    %legend_str{end+1} = 'Mars mean-equatorial radius';
    h4 = polar(theta_plot,radius_mars_atmos,'-.'); 
    set(h4,'color',cc(5,:))
    %legend_str{end+1} = 'Mars atmosphere limit';
    set(gca,'Visible','off')
    set(gcf,'color',[1 1 1])
    
% plot aerocapture trajcetories
L1 = plot(out_1_c.out.R(index1,1),out_1_c.out.R(index1,2),'color',c1);

% plot kepler (till theta = pi)
theta_kep = [out_1_c.out.okep.theta:0.02:pi-0.01 pi-0.01:0.0002:pi];
rk = out_1_c.out.okep.a * (1- out_1_c.out.okep.e^2) ./ (1 + out_1_c.out.okep.e .* cos(theta_kep));
h2 = polar(theta_kep+out_1_c.out.okep.theta_p,rk); 


% plot parking orbit
theta_park = 0:0.02:2*pi;
rp = out_p.out.a_n * (1- out_p.out.e_n^2) ./ (1 + out_p.out.e_n .* cos(theta_park));
h2 = polar(theta_park+out_1_c.out.okep.theta_p,rp); 

% plot reentry kepler
theta_e = [pi:0.02:2*pi+out_p.out.theta_entry-0.1 2*pi+out_p.out.theta_entry-0.1:0.0002:2*pi+out_p.out.theta_entry];
rkepe = out_p.out.a_e * (1- out_p.out.e_e^2) ./ (1 + out_p.out.e_e .* cos(theta_e));
h3 = polar(theta_e+out_1_c.out.okep.theta_p,rkepe); 
 
% plot entry trajcetories
L2 = plot(out_1.out.R(index2,1),out_1.out.R(index2,2),'color',c1);


if (export_figures)
    matlab2tikz('.\total_orbit.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end