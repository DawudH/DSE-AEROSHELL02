%% Thermal tool FR
% Suthes & Lucas
%% Log:
% 2-6 10:48 -- Validated correctly with 1-layer copper
clear all
close all
clc

%% Define input

%lay-up (t[mm], k[w/m/K], rho[kg/m3], cp[J/kg/K])
filename   = 'VAL3.txt';
layupin    = dlmread(filename);
% SI units
layup      = zeros(size(layupin));
layup(:,1) = layupin(:,1)./1000;
layup(:,2:5) = layupin(:,2:5);
L = layup(:,1);

%Asro input, Tinf


%Aero input, qsdot
load('heatfluxinput.mat','T','t','qmax')
qearo = qmax;
Tatm  = T; 
clear('T')

% time
%ttot = t(end);  % end time of the orbit [s]
ttot = 100;  % end time of the orbit [s]
dt   = 0.001;    % time step, chooseable
nmax = int32(ttot/dt);  % number of time steps

% change spacing of input YET TO DO!


% spaceing
L = int32(round(L*10000));


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
maxdx = double(gcd_old)/10000;
L = double(L)/10000;


Ltot = sum(L);      % length of the layup [m]
fact    = 10;           % multiplication factor of number of space steps.
dx   = maxdx/fact;     % length of a space step
imax =  int32(Ltot/dx + 1); % maximum amount of points in spaceing

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

%% Implement Cranck-Nickelson

% i is space, n is time
% space is rows, time is columns
T = zeros(imax,nmax);
T0 = 293;
T(:,1) = T0;
q = ones(1,nmax)*30000;

w = zeros(imax,1);
w(1:end-1) = v/2;
w(2:end) = w(2:end) + v/2;
C = diag(1+w) - diag(0.5*v,1) - diag(0.5*v,-1);
N = diag(1-w) + diag(0.5*v,1) + diag(0.5*v,-1);

%define matrices
A = zeros(imax,1); 
Cinv = inv(sparse(C));
G = full(Cinv*sparse(N));

for n=1:nmax-1
    A(1) = (q(n)/k(1))*v(1)*dx;
    H = full(Cinv*sparse(A));
    T(:,n+1) = G*T(:,n) + H;
end
%% Hello
t = [0:double(nmax-1)]*dt;
x = [0:double(imax-1)]*dx;
figure;
hold on
contourf(t,x,T);
colormap parula
colorbar
% title('Temperature over time and distance','Interpreter','LaTex','FontSize',13)
% xlabel('Time [s]','Interpreter','LaTex','FontSize',14)
% ylabel('Depth [m]','Interpreter','LaTex','FontSize',14)

for j = 2:length(indexx)-1
    plot([0,nmax*dt],[indexx(j)*dx,indexx(j)*dx],'--','Color', [255.0/256.0,130.0/256.0,28.0/256.0])
end
% x = L;
% y = 0:ttot:dt;
% z = ;
% contourf('x','y','z')
    
    
    
    
    
