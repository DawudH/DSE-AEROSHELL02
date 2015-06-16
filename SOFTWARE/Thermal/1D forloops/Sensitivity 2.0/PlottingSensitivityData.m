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

% Second plot data: Area(1)
%Lay2Aj = dlmread('Layup2-A-justorbit.txt');
%Lay2Ao = dlmread('Layup2-A-onetime.txt');

% Third plot data: Area(2)
%Lay4Aj = dlmread('Layup4-A-justorbit.txt');
%Lay4Ao = dlmread('Layup4-A-onetime.txt');

% Maybe 4th plot: time

%% writing in more sensible data

fluxfactor = Lay3qj(:,1);
%A          = Lay3Aj(:,1);
%timefactor

L1qj = Lay1qj(:,2);
L2qj = Lay2qj(:,2);
L3qj = Lay3qj(:,2);

%L2Aj = Lay2Aj(:,2);
%L2Ao = Lay2Ao(:,1);

%L4Aj = Lay4Aj(:,2);
%L4Ao = Lay4Ao(:,1);


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












