%% Thermal tool
% Suthes & Lucas 

clc 
clear all
close all

%% inputs
filename = 'VAL.txt';
dx = 0.005;
dt = 0.003;

%% From aero-data
% time = 0:dt:max(t);
% q = interp1(t,qmax,tempt);
% T_inf = interp1(t,T,tempt);
% totalt = max(t)/dt+1;
totalt = 40000;
% Emissivity
eps = 0.443; %[-]
time = dt*(0:totalt-1);


%% Setup input matrix (layers)
% In order of collumns [t, k, rho, cp]
% units [mm, W/m*K, kg/m3, J/kg*K]
input = dlmread(filename);
input(:,1) = input(:,1)/1000; % Layer-thickness in [m]
layers = zeros(int32(sum(input(:,1))/dx)+1,4);
l=1;
for a=1:length(input(:,1))
    bend = int32((input(a,1)/dx));
    for j=1:bend
        layers(l,1) = (l-1)*dx;
        layers(l,2:4) = input(a,2:4);
        l = l+1;
    end  
end
layers(end,:) = [layers(end-1,1)+dx,layers(end-1,2:4)];
imax = length(layers(:,1));

%% Check dt,dx
G = input(:,2)./input(:,3)./input(:,4)*dt/dx/dx;
if sum(G>0.5|G<0)>0
    disp('wrong dt or dx')
    input(:,3).*input(:,4)./input(:,2)*dx*dx/6
end
%% Input T,q:
% variable q:
x=0:totalt-1;
qmax=100000;
c1 = -25/4*qmax/totalt/totalt;
c2 = 5*qmax/totalt;
q = max(c1.*x.^2+c2.*x,0);
q = ones(1,totalt)*300000;
% variable temp:
tt = linspace(0,pi,totalt);
T_inf = 1250.*sin(tt)*0+4*0;




%% Temp-table
T0 = 20+273.15;
T = ones(imax,totalt)*(T0); %Temp-table
%T(:,1) = 300; %[K]


%% Setup FTCS-matrix
F = zeros(imax,imax);
k = layers(:,2);
rho = layers(:,3);
cp = layers(:,4);
sigma = 5.6704*10^(-8); %[W/(m^2*K^4)]
% Conducttables
Cm1 = zeros(imax,1);
Cp1 = zeros(imax,1);
Cm1(2:imax) = 2/dx./(1./k(1:imax-1)+1./k(2:imax));
Cp1(1:imax-1) = 2/dx./(1./k(2:imax)+1./k(1:imax-1));



% v = dt./rho./cp/dx;
% F = F + diag(1-v.*(Cm1+Cp1)) + diag(v(1:end-1).*Cp1(1:end-1),1) + diag(v(2:end).*Cm1(2:end),-1);
% F(1,1:2) = [1,0]+[-1,1]*v(1)*Cp1(1);
% F(imax,imax-1:imax) = [0,1]+[1,-1]*v(imax)*Cm1(imax);
alpha=k./rho./cp;
v = alpha*dt/dx/dx;
F = diag(1-2*v) + diag(v(1:end-1),1) + diag(v(2:end),-1);
F(1,1:2) = [1,0]+[-1,1]*v(1);
F(imax,imax-1:imax) = [0,1]+[1,-1]*v(imax);
%% Perform FCTS
A = zeros(imax,1);

for n=1:totalt-1
    cTrad = T(1,n)^4 - T_inf(n)^4;
    A(1,1) =  2*v(1)*dx/k(1)*(q(n)-0*eps*sigma*cTrad);
    T(:,n+1) = F*T(:,n) + A;
end


%% Plot x-t-T
plotting = 0;
if plotting
    figure;
    x = layers(:,1);
    y = time;
    z = T.';
    contourf(x,y,z)
    colormap parula
    colorbar
    title('Temperature over time and distance','Interpreter','LaTex','FontSize',13)
    xlabel('Length [m]','Interpreter','LaTex','FontSize',14)
    ylabel('Time [s]','Interpreter','LaTex','FontSize',14)

    hold on
    axis([0 layers(end,1) 0 time(end)]);
    i = length(input(:,1));
    indices = int32([1;cumsum(input(:,1)/dx)+1]);
    indices(end) = indices(end)-1;
    for i=2:i
        plot([double(indices(i)-1)*dx,double(indices(i)-1)*dx],[0,1.1*time(end)],'--','Color', [136.0/256.0,0.0,21.0/256.0]);    
    end
end
%sum(input(:,1).*input(:,3))
%alpha=k(1)/rho(1)/cp(1);tau=t*dt;k=k(1);x=dx*(imax-1);
%2*q(1)*sqrt(alpha*tau/pi)/k*exp(-x*x/4/alpha/tau)-q(1)*x/k*erfc(x/2/sqrt(alpha*tau))

%% Validation (single layer)
ssana = q(1)/k(1)*dx*(imax-1);
ssnum = T(1,end)-T(end,end);
alpha=k(1)/rho(1)/cp(1);tau=totalt*dt;k0=k(1);x=dx*(imax-1);q0=q(1);
usssana = 2*q0*sqrt(alpha*tau/pi)/k0;
usssnum = T(1,end)-T(1,1);
ussbana = 2*q0*sqrt(alpha*tau/pi)/k0*exp(-x^2/4/alpha/tau)-q0*x/k0*erfc(x/2/sqrt(alpha*tau));
ussbnum = T(end,end)-T(1,1);
output = table([ssana;ssnum;(ssnum-ssana)/ssana*100],[usssana;usssnum;(usssnum-usssana)/usssana*100],[ussbana;ussbnum;(ussbnum-ussbana)/ussbana*100],'RowNames',{'Analytical','Numerical','Difference%'},'VariableNames',{'SteadySB','UnsteadySS','UnsteadySB'})
q0/k0*sqrt(pi*alpha*tau)


qsb = 2*q0*sqrt(alpha*time/pi)/k0.*exp(-x^2/4/alpha./time)-q0*x/k0.*erfc(x/2./sqrt(alpha.*time));
x=0;
qss = 2*q0*sqrt(alpha*time/pi)/k0.*exp(-x^2/4/alpha./time)-q0*x/k0.*erfc(x/2./sqrt(alpha.*time));