%% Thermal tool FR
% Suthes & Lucas
% v1 - Implement contact resistance
%% Log:
% 2-6 10:48 -- Validated correctly with 1-layer copper
clear all
close all
clc

%% Define input

%lay-up (t[mm], k[w/m/K], rho[kg/m3], cp[J/kg/K])
filename   = 'layup3.txt';
layupin    = dlmread(filename);
% SI units
layup      = zeros(size(layupin));
layup(:,1) = layupin(:,1)./1000;
layup(:,2:6) = layupin(:,2:6);
L = layup(:,1);

% Aero input, qsdot
load('heatflux.mat','T','t','qmax_array')
tq = find(not(qmax_array==0));
qaero = qmax_array(tq(1):tq(end));
Tatm  = T(tq(1):tq(end)); 
timeq = t(tq(1):tq(end))-t(tq(1));
clear('T','t')


% time and aeroheat
%ttot = t(end);  % end time of the orbit [s]
ttot = timeq(end);  % end time of the orbit [s]
% tdur = 90;
dt   = 0.1;    % time step, chooseable
T0 = 293;
q0 = qaero*10000;
nmax = int32(ttot/dt);  % number of time steps
t = [0:double(nmax-1)]*dt;
fact    = 1;           % multiplication factor of number of space steps.
kfact = [2.5e-4;2.5e-4;2.5e-5;2.5e-5;2.5e-5];


% spaceing
L = int32(round(L*10000000));


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
maxdx = double(gcd_old)/10000000;
L = double(L)/10000000;


Ltot = sum(L);      % length of the layup [m]
dx   = maxdx/fact;     % length of a space step
imax =  int32(Ltot/dx + 1); % maximum amount of points in spaceing

%Determination of radiation parameters
emiss = layup(1,5);
sig   = 5.670373e-8; %[W/m2/K4]

%% Material properties
layers = length(L)-1;
k = zeros(imax+layers-1,1);
rho = zeros(imax+layers-1,1);
cp = zeros(imax+layers-1,1);
x = zeros(imax+layers,1);

%determine material layer thicknesses
indexx = zeros(length(L)+1,1);
indexx(1) = 1;
indexx(2:end) = int32(cumsum(L/dx))+1;
for j = 1:length(L);
    k(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,2);
    rho(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,3);
    cp(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,4);
    x(indexx(j)+j-1:indexx(j+1)+j-1) = dx*[indexx(j)-1:indexx(j+1)-1];
end
for j = 1:layers
    k(indexx(j+1)+j-1) = kfact(j)*(k(indexx(j+1)+j-2)+k(indexx(j+1)+j))/(k(indexx(j+1)+j-2)*k(indexx(j+1)+j));
    rho(indexx(j+1)+j-1) = (rho(indexx(j+1)+j-2)+rho(indexx(j+1)+j))/2;
    cp(indexx(j+1)+j-1) = (cp(indexx(j+1)+j-2)+cp(indexx(j+1)+j))/2;
end
alpha = k./rho./cp;
v = alpha*dt/dx/dx;


%% Implement Cranck-Nickelson with voids

% i is space, n is time
% space is rows, time is columns
T = zeros(imax+layers,nmax);
T(:,1) = T0;
qs = interp1(timeq,qaero,t)*10000;
Tamb = interp1(timeq,Tatm,t);


w = zeros(imax+layers,1);
w(1:end-1) = v/2;
w(2:end) = w(2:end) + v/2;
C = diag(1+w) - diag(0.5*v,1) - diag(0.5*v,-1);
N = diag(1-w) + diag(0.5*v,1) + diag(0.5*v,-1);

%define matrices
A = zeros(imax+layers,1); 
Cinv = inv(sparse(C));
G = full(Cinv*sparse(N));

for n=1:nmax-1
    qrs = -emiss*sig*(T(1,n)^4-Tamb(n)^4);
    qrb = -emiss*sig*(T(end,n)^4-Tamb(n)^4);
    A(1) = ((qs(n)+qrs)/k(1))*v(1)*dx;
    A(end) = (qrb/k(end))*v(end)*dx;
    H = full(Cinv*sparse(A));
    T(:,n+1) = G*T(:,n) + H;
end

% 
% %% Contour-Plot
% contourplot = 0;
% if contourplot
%     % t = [0:double(nmax-1)]*dt;
%     x = [0:double(imax-1)]*dx;
%     figure;
%     hold on
%     contourf(t,x,T);
%     colormap parula
%     colorbar
%     % title('Temperature over time and distance','Interpreter','LaTex','FontSize',13)
%     % xlabel('Time [s]','Interpreter','LaTex','FontSize',14)
%     % ylabel('Depth [m]','Interpreter','LaTex','FontSize',14)
% 
%     for j = 2:length(indexx)-1
%         plot([0,nmax*dt],[indexx(j)*dx,indexx(j)*dx],'--','Color', [255.0/256.0,130.0/256.0,28.0/256.0])
%     end
% end
% %% Layer analysis output
% results = zeros(length(indexx)-1,1);
% layernames = cell(length(indexx)-1,1);
% for j=1:length(indexx)-1
%     results(j) = max(T(indexx(j)+1,:));
%     layernames{j} = strcat('Layer',int2str(j));
% end
% output = table(L*1000,results,layup(:,2),layup(:,3),layup(:,4),'RowNames',layernames,'VariableNames',{'Thickness','maxT','k','rho','cp'});
%

%% Contour Plot
contourplot = 1;
if contourplot
    xcont = x;
    for j = 1:length(L)-1
        xcont(indexx(j+1)+j-1:indexx(j+1)+j) = [x(indexx(j+1)+j-1)-dx/100,x(indexx(j+1)+j)+dx/100];
    end
    figure;
    hold on
    contourf(t,xcont,T)
    colormap parula
	colorbar
    title('Temperature over time and distance')
    xlabel('Time [s]')
    ylabel('Depth [m]')
    for j = 1:length(L)-1
        plot([0,nmax*dt],[x(indexx(j+1)+j),x(indexx(j+1)+j)],'--','Color', [255.0/256.0,130.0/256.0,28.0/256.0])
    end
    axis ij   
end


%% Layer analysis
results = zeros(length(L),1);
layernames = cell(length(L),1);
for j=1:length(L)
    results(j) = max(max(T(indexx(j)+j-1:indexx(j+1)+j-1,:)));
    layernames{j} = strcat('Layer',int2str(j));
end
output = table(layup(:,1)*1000,results,layup(:,6),layup(:,2),layup(:,3),layup(:,4),'RowNames',layernames,'VariableNames',{'Thickness','maxT','allowT','k','rho','cp'});


%% Validation
valid = 0;  
if valid
    valres = dlmread('layup1res.txt');
    plotS = zeros(nmax,8);
    plotS(:,1) = T(1,:).';
    for j = 2:length(indexx)-1
        plotS(:,2*j-2:2*j-1) = [T(indexx(j)+j-3,:).',T(indexx(j)+j-1,:).'];
    end
    plotS(:,8) = T(end,:).';
    plotVAL = zeros(nmax,8);
    for j = 2:length(valres(1,:))
        plotVAL(:,j-1) = interp1(valres(:,1),valres(:,j),t).';
    end
%     figure;
%     hold on
%     plot(t,plotS,'--')
%     ax = gca;
%     ax.ColorOrderIndex = 1;
%     plot(t,plotVAL)
%     figure;
%     plot(t,abs(plotS-plotVAL)./plotVAL)
%     max(abs(plotS-plotVAL)./plotVAL)
end
    
disp(output)    
    
    
    
