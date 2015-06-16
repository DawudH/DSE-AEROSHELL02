%% Plotting al the sensitivity data
% Suthes & Lucas
clear all
close all
clc

% For plotting in LaTeX style
addpath('..\..\..\matlab2tikz')

%% Reading in the sensitivity data

% First plot data: flux
Lay2qj = dlmread('Layup2-q_q0-justorbit.txt');
Lay4qj = dlmread('Layup4-q_q0-justorbit.txt');

% Second plot data: Area(1)
Lay2Aj = dlmread('Layup2-A-justorbit.txt');
Lay2Ao = dlmread('Layup2-A-onetime.txt');

% Third plot data: Area(2)
Lay4Aj = dlmread('Layup4-A-justorbit.txt');
Lay4Ao = dlmread('Layup4-A-onetime.txt');

% Maybe 4th plot: time

%% writing in more sensible data

fluxfactor = Lay2qj(:,1);
A          = Lay2Aj(:,1);
%timefactor

L2qj = Lay2qj(:,2);
L4qj = Lay4qj(:,2);

L2Aj = Lay2Aj(:,2);
L2Ao = Lay2Ao(:,1);

L4Aj = Lay4Aj(:,2);
L4Ao = Lay4Ao(:,1);


%% Plotting the obtained data in LaTeX style


cc = parula(4);
figure(1)
plot(fluxfactor,L2qj,'d-','color',cc(1,:))
hold on
plot(fluxfactor,L4qj,'+-','color',cc(2,:))
grid on
xlabel('Heat flux ratio [-]')
ylabel('Areal mass [kg/m^2 ]')
legend('Layup 2','Layup 4','Location','northwest')
matlab2tikz('.\LaTeX\fluxsensitivity.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);












