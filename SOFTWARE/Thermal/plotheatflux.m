clear all
close all
clc

% For plotting in LaTeX style
addpath('..\matlab2tikz')



% Aero input, qsdot
load('heatflux_input/aerocapture_rho_1.mat','T','t','qmax_array','A_whetted')
tq1 = find(not(qmax_array<0.01));
qaero1 = qmax_array(tq1(1):tq1(end));
Tatm1  = T(tq1(1):tq1(end)); 
timeq1 = t(tq1(1):tq1(end))-t(tq1(1));
clear('T','t')
% Aero input, qsdot
load('heatflux_input/aerocapture_rho_0_9.mat','T','t','qmax_array','A_whetted')
tq2 = find(not(qmax_array<0.01));
qaero2 = qmax_array(tq2(1):tq2(end));
Tatm2  = T(tq2(1):tq2(end)); 
timeq2 = t(tq2(1):tq2(end))-t(tq2(1));
clear('T','t')
% Aero input, qsdot
load('heatflux_input/aerocapture_rho_1_1.mat','T','t','qmax_array','A_whetted')
tq3 = find(not(qmax_array<0.01));
qaero3 = qmax_array(tq3(1):tq3(end));
Tatm3  = T(tq3(1):tq3(end)); 
timeq3 = t(tq3(1):tq3(end))-t(tq3(1));
clear('T','t')

d = 200;
e = d/4;
figure
hold on
cc = parula(5);
plot(timeq2(1),qaero2(1),'+-','color',cc(2,:))
plot(timeq1(1),qaero1(1),'d-','color',cc(1,:))
plot(timeq3(1),qaero3(1),'s-','color',cc(3,:))
plot(timeq2(1:d:end),qaero2(1:d:end),'+','color',cc(2,:))
plot(timeq1(1:d:end),qaero1(1:d:end),'d','color',cc(1,:))
plot(timeq3(1:d:end),qaero3(1:d:end),'s','color',cc(3,:))
plot(timeq2(1:e:end),qaero2(1:e:end),'-','color',cc(2,:))
plot(timeq1(1:e:end),qaero1(1:e:end),'-','color',cc(1,:))
plot(timeq3(1:e:end),qaero3(1:e:end),'-','color',cc(3,:))
grid on
ylabel('$\dot{q}_s$ $\left[ W\cdot cm^{-2} \right]$','Interpreter','LaTeX')
xlabel('t $\left[ s \right]$','interpreter','LaTeX')
legend('\rho/\rho_0 = 0.9','\rho/\rho_0 = 1.0','\rho/\rho_0 = 1.1','Location','northeast')
matlab2tikz('.\Figures\plotheatflux.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

