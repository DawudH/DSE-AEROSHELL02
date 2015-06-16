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

%% processing (plot/write to file)
figure('name','parameters over time')
t_1 = out_1.out.tp(1:end-1);
t_1_1 = out_1_1.out.tp(1:end-1);
t_0_9 = out_0_9.out.tp(1:end-1);

cc = parula(7);
c1 = cc(1,:);
c2 = cc(3,:);
c3 = cc(5,:);

xlimits = [0 max(t_0_9)];

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
    matlab2tikz('.\orbit_sensitivity.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
end