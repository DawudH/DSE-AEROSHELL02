%% Thermal tool FR
% Suthes & Lucas
% v1 - Successfully validated using copper-example from paper
% v2 - Implemented FTCS-matrix
% v3 - Implemented BTCS-matrix
% v4 - Implemented Crank-Nicolson

clc 
clear all
close all

% For plotting in LaTeX style
addpath('..\matlab2tikz')


nmax = 10000;
imax = 1001;%10000
k = ones(imax,1)*386.0;
rho = ones(imax,1)*8954.0;
cp = ones(imax,1)*383.1;

dx = 0.0005;%0.00005
dt = 0.1;

alpha = k./rho./cp;
v = alpha*dt/dx/dx;

% i is space, n is time
% space is rows, time is columns
%T = zeros(imax,1);
%U = zeros(2,nmax);
T = zeros(imax,nmax);
T0 = 293;
T(:,1) = T0;
%U(:,1) = T0;
q = ones(1,nmax)*300000;

C = diag(1+v) - 0.5*diag(v(1:end-1),1) - 0.5*diag(v(2:end),-1);
C(1,1) = 1+0.5*v(1);
C(imax,imax) = 1+0.5*v(imax);
N = diag(1-v) + 0.5*diag(v(1:end-1),1) + 0.5*diag(v(2:end),-1);
N(1,1) = 1-0.5*v(1);
N(imax,imax) = 1-0.5*v(imax);

A = zeros(imax,1);
invC = inv(C);
G = invC*N;
for n=1:nmax-1
    A(1) = (q(n)/k(1))*v(1)*dx;
    H = invC*A;
    T(:,n+1) = G*T(:,n) + H;
    %T = G*T + H;
    %U(:,n+1) = [T(1),T(end)];
end
size(T)
q = q(1);
k = k(1);
alpha = alpha(1);

xlist = (0:imax/20:imax-1)*dx;
tlist = (0:nmax-1)*dt;
xlist = 0;
valTs = T0 + (2*q.*sqrt(alpha.*tlist/pi)/k).*exp(-xlist.*xlist./(4*alpha.*tlist)) - (q.*xlist/k).*(1-erf(xlist./(2.*sqrt(alpha.*tlist))));
xlist = 0.25;
valTm = T0 + (2*q.*sqrt(alpha.*tlist/pi)/k).*exp(-xlist.*xlist./(4*alpha.*tlist)) - (q.*xlist/k).*(1-erf(xlist./(2.*sqrt(alpha.*tlist))));
xlist = 0.5;
valTb = T0 + (2*q.*sqrt(alpha.*tlist/pi)/k).*exp(-xlist.*xlist./(4*alpha.*tlist)) - (q.*xlist/k).*(1-erf(xlist./(2.*sqrt(alpha.*tlist))));
figure;
hold on
cc = parula(5);
plot(tlist(1),T(1,1),'+--','color',cc(1,:))
plot(tlist(1),valTs(1),'+-','color',cc(1,:))
plot(tlist(1),T(501,1),'d--','color',cc(2,:))
plot(tlist(1),valTm(1),'d-','color',cc(2,:))
plot(tlist(1),T(end,1),'s--','color',cc(3,:))
plot(tlist(1),valTb(1),'s-','color',cc(3,:))
plot(tlist(1:500:end),T(1,1:500:end),'+','color',cc(1,:))
plot(tlist(1:500:end),valTs(1:500:end),'+','color',cc(1,:))
plot(tlist(1:500:end),T(501,1:500:end),'d','color',cc(2,:))
plot(tlist(1:500:end),valTm(1:500:end),'d','color',cc(2,:))
plot(tlist(1:500:end),T(end,1:500:end),'s','color',cc(3,:))
plot(tlist(1:500:end),valTb(1:500:end),'s','color',cc(3,:))
plot(tlist(1:50:end),T(1,1:50:end),'--','color',cc(1,:))
plot(tlist(1:50:end),valTs(1:50:end),'-','color',cc(1,:))
plot(tlist(1:50:end),T(501,1:50:end),'--','color',cc(2,:))
plot(tlist(1:50:end),valTm(1:50:end),'-','color',cc(2,:))
plot(tlist(1:50:end),T(end,1:50:end),'--','color',cc(3,:))
plot(tlist(1:50:end),valTb(1:50:end),'-','color',cc(3,:))
grid on
ylabel('T $\left[ K \right]$','Interpreter','LaTeX')
xlabel('t $\left[ s \right]$','interpreter','LaTeX')
legend('x = 0.0, Numerical','x = 0.0, Analytical','x = 0.5, Numerical','x = 0.5, Analytical','x = 1.0, Numerical','x = 1.0, Analytcal','Location','northwest')

% x = 0;
% t = nmax*dt;
% anas = T0 + (2*q*sqrt(alpha*t/pi)/k)*exp(-x*x/(4*alpha*t)) - (q*x/k)*(1-erf(x/(2*sqrt(alpha*t))));
% x = (imax-1)*dx;
% anab = T0 + (2*q*sqrt(alpha*t/pi)/k)*exp(-x*x/(4*alpha*t)) - (q*x/k)*(1-erf(x/(2*sqrt(alpha*t))));
% 
% output = table([T0;T(1,1);(T(1,1)-T0)/T0*100],[T0;T(imax,1);(T(imax,1)-T0)/T0*100],[anas;T(1,nmax);(T(1,nmax)-anas)/anas*100],[anab;T(imax,nmax);(T(imax,nmax)-anab)/anab*100],'RowNames',{'Analytical','Numerical','Difference'},'VariableNames',{'Surface','Back','Surfaceend','Backend'});
% 
% output
% 
% xlist = 1:imax/20:imax;
% numlist120 = T(xlist,nmax);
% realx = (xlist-1)*dx;
% analist120 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% t = nmax*0.8*dt;
% numlist96 = T(xlist,nmax*0.8);
% analist96 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% t = nmax*0.6*dt;
% numlist72 = T(xlist,nmax*0.6);
% analist72 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% t = nmax*0.4*dt;
% numlist48 = T(xlist,nmax*0.4);
% analist48 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% t = nmax*0.2*dt;
% numlist24 = T(xlist,nmax*0.2);
% analist24 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% t = 0;
% numlist0 = T(xlist,1);
% analist0 = T0 + (2*q*sqrt(alpha*t/pi)/k).*exp(-realx.*realx/(4*alpha*t)) - (q.*realx/k).*(1-erf(realx/(2*sqrt(alpha*t))));
% 
% figure;
% hold on
% plot(analist120,realx,'b-',numlist120,realx,'r^')
% plot(analist96,realx,'b-',numlist96,realx,'r^')
% plot(analist72,realx,'b-',numlist72,realx,'r^')
% plot(analist48,realx,'b-',numlist48,realx,'r^')
% plot(analist24,realx,'b-',numlist24,realx,'r^')
% plot(analist0,realx,'b-',numlist0,realx,'r^')
% axis ij

