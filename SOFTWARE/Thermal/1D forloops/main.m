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
load('heatfluxinput.mat')
qearo = qmax;
Tatm  = T;

% time
ttot = t(end);  % end time of the orbit [s]
dt   = 0.01;    % time step, chooseable
nmax = (ttot/dt);  % number of time steps


% spaceing
L = L*10000;

if length(L)==1
    gcd_old = L(1);
else
    gcd_old = gcd(L(1),L(2));
    if length(L)>2
        for  j = 3:length(L)
            gcd_old = gcd(gcd_old,L(j));
        end
    end
end
maxdx = gcd_old/10000;
L = L/10000;


Ltot = sum(L);      % length of the layup [m]
N    = 10;           % multiplication factor of number of space steps.
dx   = maxdx/N;     % length of a space step
imax =  Ltot/dx +1; % maximum amount of points in spaceing

%Determination of radiation parameters
emiss = layup(1,5);
sig   = 5.670373e-8; %[W/m2/K4]

%% Material properties
k     = zeros(imax-1,1); 
rho   = zeros(imax-1,1);
cp    = zeros(imax-1,1);

%determine material layer thicknesses
indexx = zeros(length(L)+1,1);
indexx(1) = 0;
indexx(2:end) = int32(cumsum(L/dx));
for j = 1:length(L);
    k(indexx(j)+1:indexx(j+1)) = layup(j,2);
    rho(indexx(j)+1:indexx(j+1)) = layup(j,3);
    cp(indexx(j)+1:indexx(j+1)) = layup(j,4);
end
alpha = k./rho./cp;
v = alpha*dt/dx/dx;



    
  
    
    
    
    
    
