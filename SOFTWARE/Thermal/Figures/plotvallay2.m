%% Thermal tool FR
% Suthes & Lucas
%% Log:
% 2-6 10:48 -- Validated correctly with 1-layer copper
clear all
close all
clc

%% Define input

%lay-up (t[mm], k[w/m/K], rho[kg/m3], cp[J/kg/K])
filename   = 'layup2.txt';
layupin    = dlmread(filename);
% SI units
layup      = zeros(size(layupin));
layup(:,1) = layupin(:,1)./1000;
layup(:,2:5) = layupin(:,2:5);
L = layup(:,1);


% time and aeroheat
%ttot = t(end);  % end time of the orbit [s]
ttot = 200;  % end time of the orbit [s]
tdur = 90;
dt   = 0.1;    % time step, chooseable
T0 = 293;
q0 = 62000;
nmax = int32(ttot/dt);  % number of time steps
t = [0:double(nmax-1)]*dt;
fact    = 1;           % multiplication factor of number of space steps.
%kfact = [2.5e-5/fact;2.5e-6/fact;2.5e-2/fact];
kfact = [2.5e-5;1;1];

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
voidlayers = length(L)-1;
kvoid = zeros(imax+voidlayers-1,1);
rhovoid = zeros(imax+voidlayers-1,1);
cpvoid = zeros(imax+voidlayers-1,1);
k     = zeros(imax-1,1); 
rho   = zeros(imax-1,1);
cp    = zeros(imax-1,1);

%determine material layer thicknesses
indexx = zeros(length(L)+1,1);
indexx(1) = 1;
indexx(2:end) = int32(cumsum(L/dx))+1;
for j = 1:length(L);
    kvoid(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,2);
    rhovoid(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,3);
    cpvoid(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,4);
    k(indexx(j):indexx(j+1)-1) = layup(j,2);
    rho(indexx(j):indexx(j+1)-1) = layup(j,3);
    cp(indexx(j):indexx(j+1)-1) = layup(j,4);
end
for j = 1:voidlayers
    kvoid(indexx(j+1)+j-1) = kfact(j)*(kvoid(indexx(j+1)+j-2)+kvoid(indexx(j+1)+j))/(kvoid(indexx(j+1)+j-2)*kvoid(indexx(j+1)+j));
    rhovoid(indexx(j+1)+j-1) = (rhovoid(indexx(j+1)+j-2)+rhovoid(indexx(j+1)+j))/2;
    cpvoid(indexx(j+1)+j-1) = (cpvoid(indexx(j+1)+j-2)+cpvoid(indexx(j+1)+j))/2;
end
alphavoid = kvoid./rhovoid./cpvoid;
vvoid = alphavoid*dt/dx/dx;

alpha = k./rho./cp;
v = alpha*dt/dx/dx;



%% Implement Cranck-Nickelson

% i is space, n is time
% space is rows, time is columns
T = zeros(imax,nmax);
T(:,1) = T0;
Tamb = ones(nmax,1)*T0;
qs = ones(nmax,1)*q0;
qs(tdur/dt:end) = 0;


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
    qr = -emiss*sig*(T(1,n)^4-Tamb(n)^4);
    A(1) = ((qs(n)+qr)/k(1))*v(1)*dx;
    H = full(Cinv*sparse(A));
    T(:,n+1) = G*T(:,n) + H;
end

%% Implement Cranck-Nickelson with voids

% i is space, n is time
% space is rows, time is columns
S = zeros(imax+voidlayers,nmax);
S(:,1) = T0;
Tamb = ones(nmax,1)*T0;
qs = ones(nmax,1)*q0;
qs(tdur/dt:end) = 0;


wvoid = zeros(imax+voidlayers,1);
wvoid(1:end-1) = vvoid/2;
wvoid(2:end) = wvoid(2:end) + vvoid/2;
C = diag(1+wvoid) - diag(0.5*vvoid,1) - diag(0.5*vvoid,-1);
N = diag(1-wvoid) + diag(0.5*vvoid,1) + diag(0.5*vvoid,-1);

%define matrices
A = zeros(imax+voidlayers,1); 
Cinv = inv(sparse(C));
G = full(Cinv*sparse(N));

for n=1:nmax-1
    qrs = -emiss*sig*(S(1,n)^4-Tamb(n)^4);
    qrb = -emiss*sig*(S(end,n)^4-Tamb(n)^4);
    A(1) = ((qs(n)+qrs)/kvoid(1))*vvoid(1)*dx;
    A(end) = (qrb/kvoid(end))*vvoid(end)*dx;
    H = full(Cinv*sparse(A));
    S(:,n+1) = G*S(:,n) + H;
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
%% Validation
valid = 1;
% if valid
%     valres = dlmread('layup1res.txt');
%     figure;
%     hold on
%     plot(valres(:,1),valres(:,2),'b-')
%     plot(valres(:,1),valres(:,4),'r-')
% %     for j = 2:length(valres(1,:))
% %         plot(valres(:,1),valres(:,j))
% %     end
%     plot(t,T(1,:),'b--')
%     plot(t,T(3,:),'r--')
% end
if valid
    valres = dlmread('layup2res.txt');

    plotS = zeros(nmax,6);
    for j = 2:length(indexx)-1
        plotS(:,2*j-3:2*j-2) = [S(indexx(j)+j-3,:).',S(indexx(j)+j-1,:).'];
    end
    plotVAL = zeros(nmax,6);
    for j = 5:length(valres(1,:))
        plotVAL(:,j-4) = interp1(valres(1:end-1,1),valres(1:end-1,j),t).';
    end
    plotVALex = zeros(nmax,3);
    for j = 2:4
        plotVALex(:,j-1) = interp1(valres(1:end-1,1),valres(1:end-1,j),t).';
    end
%     figure;
%     hold on   
%     
%     plot(t,plotS,'--')
%     ax = gca;
%     ax.ColorOrderIndex = 1;
%     plot(t,plotVAL)
%     figure;
%     plot(t,abs(plotS-plotVAL)./plotVAL)
%     max(abs(plotS-plotVAL)./plotVAL)
end
%o,x,*,s,d,+
figure;


subplot(1,2,1)
hold on
cc = parula(4);

plot(t(1),plotVALex(1,1),'+-','color',cc(1,:))
plot(t(1),plotVALex(1,2),'d-','color',cc(2,:))
plot(t(1),plotVALex(1,3),'s-','color',cc(3,:))
plot(t(1),(plotS(1,1)+plotS(1,2))/2,'+--','color',cc(1,:))
plot(t(1),(plotS(1,3)+plotS(1,4))/2,'d--','color',cc(2,:))
plot(t(1),(plotS(1,5)+plotS(1,6))/2,'s--','color',cc(3,:))
d = 200;
e = 40;
plot(t(1:e:end),plotVALex(1:e:end,1),'-','color',cc(1,:))
plot(t(1:e:end),(plotS(1:e:end,1)+plotS(1:e:end,2))/2,'--','color',cc(1,:))
plot(t(1:e:end),plotVALex(1:e:end,2),'-','color',cc(2,:))
plot(t(1:e:end),(plotS(1:e:end,3)+plotS(1:e:end,4))/2,'--','color',cc(2,:))
plot(t(1:e:end),plotVALex(1:e:end,3),'-','color',cc(3,:))
plot(t(1:e:end),(plotS(1:e:end,5)+plotS(1:e:end,6))/2,'--','color',cc(3,:))

plot(t(1:d:end),plotVALex(1:d:end,1),'+','color',cc(1,:))
plot(t(1:d:end),(plotS(1:d:end,1)+plotS(1:d:end,2))/2,'+','color',cc(1,:))
plot(t(1:d:end),plotVALex(1:d:end,2),'d','color',cc(2,:))
plot(t(1:d:end),(plotS(1:d:end,3)+plotS(1:d:end,4))/2,'d','color',cc(2,:))
plot(t(1:d:end),plotVALex(1:d:end,3),'s','color',cc(3,:))
plot(t(1:d:end),(plotS(1:d:end,5)+plotS(1:d:end,6))/2,'s','color',cc(3,:))
ylabel('T $\left[ K \right]$','Interpreter','LaTeX')
xlabel('t $\left[ s \right]$','interpreter','LaTeX')
legend('TC1 Experimental','TC2 Experimental','TC3 Experimental','TC1 Model','TC2 Model','TC3 Model','Location','northeast')



subplot(1,2,2)
hold on
plot(t(1),100*abs((plotS(1,1)+plotS(1,2))/2-plotVALex(1,1))./plotVALex(1,1),'+-','color',cc(1,:))
plot(t(1),100*abs((plotS(1,3)+plotS(1,4))/2-plotVALex(1,2))./plotVALex(1,2),'d-','color',cc(2,:))
plot(t(1),100*abs((plotS(1,5)+plotS(1,6))/2-plotVALex(1,3))./plotVALex(1,3),'s-','color',cc(3,:))


plot(t(1:e:end),100*abs((plotS(1:e:end,1)+plotS(1:e:end,2))/2-plotVALex(1:e:end,1))./plotVALex(1:e:end,1),'-','color',cc(1,:))
plot(t(1:e:end),100*abs((plotS(1:e:end,3)+plotS(1:e:end,4))/2-plotVALex(1:e:end,2))./plotVALex(1:e:end,2),'-','color',cc(2,:))
plot(t(1:e:end),100*abs((plotS(1:e:end,5)+plotS(1:e:end,6))/2-plotVALex(1:e:end,3))./plotVALex(1:e:end,3),'-','color',cc(3,:))


plot(t(1:d:end),100*abs((plotS(1:d:end,1)+plotS(1:d:end,2))/2-plotVALex(1:d:end,1))./plotVALex(1:d:end,1),'+','color',cc(1,:))
plot(t(1:d:end),100*abs((plotS(1:d:end,3)+plotS(1:d:end,4))/2-plotVALex(1:d:end,2))./plotVALex(1:d:end,2),'d','color',cc(2,:))
plot(t(1:d:end),100*abs((plotS(1:d:end,5)+plotS(1:d:end,6))/2-plotVALex(1:d:end,3))./plotVALex(1:d:end,3),'s','color',cc(3,:))

ylabel('Error $\left[ \% \right]$','Interpreter','LaTeX')
xlabel('t $\left[ s \right]$','interpreter','LaTeX')
%legend('TC1','TC2','TC3','Location','northwest')    
addpath('..\matlab2tikz')

matlab2tikz('.\Figures\plotvallay2.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
    
    
