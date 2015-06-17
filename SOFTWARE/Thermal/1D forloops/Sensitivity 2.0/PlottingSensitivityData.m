%% Plotting al the sensitivity data
% Suthes & Lucas
clear all
close all
clc

% For plotting in LaTeX style
addpath('..\..\..\matlab2tikz')

%% Reading in the sensitivity data

% First plot data: flux
Lay1qj = dlmread('Layup1-q_q0-justorbit.txt');
Lay2qj = dlmread('Layup2-q_q0-justorbit.txt');
Lay3qj = dlmread('Layup3-q_q0-justorbit.txt');

% Area orbit
Lay1Aj = dlmread('Layup1-A-justorbit.txt');
Lay2Aj = dlmread('Layup2-A-justorbit.txt');
Lay3Aj = dlmread('Layup3-A-justorbit.txt');

% Area onetime
Lay1Ao = dlmread('Layup1-A-onetime.txt');
Lay2Ao = dlmread('Layup2-A-onetime.txt');
Lay3Ao = dlmread('Layup3-A-onetime.txt');

% Maybe 4th plot: time



%% Plotting the obtained data in LaTeX style


cc = parula(5);
figure(1)
plot(Lay1qj(:,1),Lay1qj(:,2),'d-','color',cc(1,:))
hold on
plot(Lay2qj(:,1),Lay2qj(:,2),'+-','color',cc(2,:))
plot(Lay3qj(:,1),Lay3qj(:,2),'s-','color',cc(3,:))
grid on
xlabel('Heat flux ratio $\left[-\right]$','Interpreter','LaTeX')
ylabel('Area density $\left[ kg \cdot m^{-2}\right]$','interpreter','latex')
legend('Layup 1','Layup 2','Layup 3','Location','northeast')
matlab2tikz('.\LaTeX\fluxsensitivity.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

cc = parula(4);
figure(2)
hold on
plot(Lay1Aj(:,1),Lay1Aj(:,3),'d-','color',cc(1,:))
plot(Lay1Ao(:,1),Lay1Ao(:,3),'d--','color',cc(1,:))
plot(Lay2Aj(:,1),Lay2Aj(:,3),'+-','color',cc(2,:))
plot(Lay2Ao(:,1),Lay2Ao(:,3),'+--','color',cc(2,:))
plot(Lay3Aj(:,1),Lay3Aj(:,3),'s-','color',cc(3,:))
plot(Lay3Ao(:,1),Lay3Ao(:,3),'s--','color',cc(3,:))
grid on
xlabel('Diameter $\left[ m \right]$','Interpreter','LaTeX')
ylabel('Frontal TPS mass $\left[ kg \right]$','interpreter','latex')
legend('Layup 1 orbit ','Layup 1 direct','Layup 2 orbit','Layup 2 direct','Layup 3 orbit','Layup 3 direct','Location','southeast')
axis([0 26 0 650])
matlab2tikz('.\LaTeX\Areasensitivity.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);












