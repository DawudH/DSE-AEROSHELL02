%% Thermal tool FR
% Suthes & Lucas 

clc 
clear all
close all

nmax = 10000;
imax = 100;
k = 386.0;
rho = 8954.0;
cp = 383.1;

dx = 0.005;
dt = 0.012;
%dt = 0.0370280375647668487;

alpha = k/rho/cp;
v = alpha*dt/dx/dx;

T = zeros(imax,nmax);
T0 = 293;
T(:,1) = 293;
q = 300000;

for n=1:nmax-1
    fac = 1;
    T(1,n+1) = (1-fac*v)*T(1,n) + fac*v*T(2,n) + (q/k)*v*dx*fac/2;%!!!! +q/k*dx is questionable
    for i=2:imax-1
        T(i,n+1) = v*T(i-1,n) + (1-2*v)*T(i,n) + v*T(i+1,n);
    % i=imax
    T(imax,n+1) = v*T(imax-1,n) + (1-v)*T(imax,n);
    end
end

num = T(imax,nmax)-T(1,1);
x = imax*dx;
t = nmax*dt;
ana = (2*q*sqrt(alpha*t/pi)/k)*exp(-x*x/(4*alpha*t)) - (q*x/k)*(1-erf(x/(2*sqrt(alpha*t))));

output = table([T0;T(1,1);(T(1,1)-T0)/T0*100],[T0;T(imax,1);(T(imax,1)-T0)/T0*100],[],[],'RowNames',{'Analytical','Numerical','Difference'},'VariableNames',{'Surface t=0','Back t=0','Surface t=end','Back t=end'})
