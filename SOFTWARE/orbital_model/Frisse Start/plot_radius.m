clc
clear
close all

variables

addpath('..\..\Aero model\orbits')

load('out_d6_just_orbit','out')
outj6 = out;
hj6 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d9_just_orbit','out')
outj9 = out;
hj9 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d12_just_orbit','out')
outj12 = out;
hj12 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d15_just_orbit','out')
outj15 = out;
hj15 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d18_just_orbit','out')
outj18 = out;
hj18 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;

load('out_d6_one_time','out')
outo6 = out;
ho6 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d9_one_time','out')
outo9 = out;
ho9 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d12_one_time','out')
outo12 = out;
ho12 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d15_one_time','out')
outo15 = out;
ho15 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;
load('out_d18_one_time','out')
outo18 = out;
ho18 = sqrt(out.R(:,1).^2 + out.R(:,3).^2 + out.R(:,2).^2) - R_m;

% figure('name','radius')
% hold on
% just=fit([6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]','cubicinterp');
% plot(just,[6,9,12,15,18],[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)])
% one=fit([6,9,12,15,18]',[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)]','cubicinterp');
% plot(one,[6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)])
% 
% figure('name','radius')
% hold on
% just = polyfit([6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]',3);
% justx = linspace(6,18,100);
% justy = polyval(just,justx);
% plot(justx,justy,[6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]','x');
% one = polyfit([6,9,12,15,18]',[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)]',3);
% onex = linspace(6,18,100);
% oney = polyval(one,onex);
% plot(onex,oney,[6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)],'d')
cc = parula(7);

figure('name','radius')

subplot(1,3,1)
hold on
grid on
xlabel('$r$ $\left[m\right]$','interpreter','latex')
ylabel('$q$ $\left[Pa\right]$','interpreter','latex')
plot([6,9,12,15,18],[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)],'-o','color',cc(1,:),'MarkerEdgeColor',cc(1,:))
plot([6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)],'-d','color',cc(3,:),'MarkerEdgeColor',cc(3,:))

% subplot(1,4,2)
% hold on
% grid on
% xlabel('$r$ $\left[m\right]$','interpreter','latex')
% ylabel('$\rho$ $\left[kg\cdot m^{-3}\right]$','interpreter','latex')
% semilogy([6,9,12,15,18],[max(outj6.rho),max(outj9.rho),max(outj12.rho),max(outj15.rho),max(outj18.rho)],'-o','color',cc(1,:),'MarkerEdgeColor',cc(1,:))
% semilogy([6,9,12,15,18],[max(outo6.rho),max(outo9.rho),max(outo12.rho),max(outo15.rho),max(outo18.rho)],'-d','color',cc(3,:),'MarkerEdgeColor',cc(3,:))

subplot(1,3,3)
hold on
grid on
xlabel('$r$ $\left[m\right]$','interpreter','latex')
ylabel('$M$ $\left[-\right]$','interpreter','latex')
ylim([40,45])
plot([6,9,12,15,18],[max(outj6.M),max(outj9.M),max(outj12.M),max(outj15.M),max(outj18.M)],'-o','color',cc(1,:),'MarkerEdgeColor',cc(1,:))
plot([6,9,12,15,18],[max(outo6.M),max(outo9.M),max(outo12.M),max(outo15.M),max(outo18.M)],'-d','color',cc(3,:),'MarkerEdgeColor',cc(3,:))

subplot(1,3,2)
hold on
grid on
xlabel('$r$ $\left[m\right]$','interpreter','latex')
ylabel('$h$ $\left[km\right]$','interpreter','latex')
%ylim([399.5,400.5])
plot([6,9,12,15,18],[min(hj6),min(hj9),min(hj12),min(hj15),min(hj18)]/1000,'-o','color',cc(1,:),'MarkerEdgeColor',cc(1,:))
plot([6,9,12,15,18],[min(ho6),min(ho9),min(ho12),min(ho15),min(ho18)]/1000,'-d','color',cc(3,:),'MarkerEdgeColor',cc(3,:))

matlab2tikz('.\LaTeX\radius.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);