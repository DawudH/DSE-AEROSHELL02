clear all
close all
clc

%% Define input

%lay-up (t[mm], k[w/m/K], rho[kg/m3], cp[J/kg/K])
filename   = 'IRVE4.txt';
layupin    = dlmread(filename);
% SI units
layup      = zeros(size(layupin));
layup(:,1) = layupin(:,1)./1000;
layup(:,2:5) = layupin(:,2:5);
L = layup(:,1);

%Asro input, Tinf


%Aero input, qsdot


% time
tmax = 120;         % end time of the orbit [s]
nt = 1000000;       % number of time steps
dt   = tmax/(nt-1);      
% space
Ltot = sum(L);      % length of the layup [m]
nx = 1200;           % number of nodes in the x direction
dx   = Ltot/(nx-1);


%Determination of radiation parameters
emiss = layup(1,5);
sig   = 5.670373e-8; %[W/m2/K4]

%% Material properties
k     = zeros(1,nx); 
rho   = zeros(1,nx);
cp    = zeros(1,nx);
alpha = zeros(1,nx);

% determine material layer thicknesses
% L = zeros(1,length(layup));
% for i = 1:length(layup);
%     L(i) = 
% end

    
    
    
    
    
    
    
