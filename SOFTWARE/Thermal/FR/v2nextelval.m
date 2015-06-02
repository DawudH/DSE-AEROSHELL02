%% Thermal tool FR
% Suthes & Lucas
% v1 - Crank-Nicolson with nextel

clc 
clear all
close all

nlevels = [10;25;40;70;100;250;400;700;1000;2500;4000;7000;10000];
results = zeros(length(nlevels),4);

for p=1:length(nlevels)
    

    nmax = nlevels(p);
    imax = 800;%10000
    k = ones(imax,1)*0.146;
    rho = ones(imax,1)*1362;
    cp = ones(imax,1)*1130.0;

    dx = 0.000025;%0.00005
    dt = 1;

    alpha = k./rho./cp;
    v = alpha*dt/dx/dx;

    % i is space, n is time
    % space is rows, time is columns
    T = zeros(imax,nmax);
    T0 = 293;
    T(:,1) = T0;
    q = ones(1,nmax)*30000;

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
    %     T(1,n+1) = (1-v)*T(1,n) + v*T(2,n) + (q/k)*v*dx;%!!!! +q/k*dx is questionable
    %     for i=2:imax-1
    %         T(i,n+1) = v*T(i-1,n) + (1-2*v)*T(i,n) + v*T(i+1,n);
    %     % i=imax
    %     T(imax,n+1) = v*T(imax-1,n) + (1-v)*T(imax,n);
    %     end
    end

    q = q(1);
    k = k(1);
    alpha = alpha(1);

    x = 0;
    t = nmax*dt;
    anas = T0 + (2*q*sqrt(alpha*t/pi)/k)*exp(-x*x/(4*alpha*t)) - (q*x/k)*(1-erf(x/(2*sqrt(alpha*t))));
    x = (imax-1)*dx;
    anab = T0 + (2*q*sqrt(alpha*t/pi)/k)*exp(-x*x/(4*alpha*t)) - (q*x/k)*(1-erf(x/(2*sqrt(alpha*t))));

    output = table([T0;T(1,1);(T(1,1)-T0)/T0*100],[T0;T(imax,1);(T(imax,1)-T0)/T0*100],[anas;T(1,nmax);(T(1,nmax)-anas)/anas*100],[anab;T(imax,nmax);(T(imax,nmax)-anab)/anab*100],'RowNames',{'Analytical','Numerical','Difference'},'VariableNames',{'Surface','Back','Surfaceend','Backend'});
    results(p,:) = [1,anas,T(1,nmax),(T(1,nmax)-anas)/anas];
end
