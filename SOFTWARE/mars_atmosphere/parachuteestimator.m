addpath('../matlab2tikz');
% load('TILL_GROUND_DO_NOT_USE_FOR_OTHER_PURPOSES_THAN_TD.mat');
% R_mars = 3.389945945211271e6; %[m]
% h = sqrt(sum(out.R.^2,2))-R_mars;
% n = out.Vm.^2./(2*h*g0);
% plot(h,n);
% axis([0, 15000, -30, 30]);

% close all;
% clear; clc;
aer = marsatmosphere();

narray = 3;%0.7:0.1:4;
totalmassarray = zeros(size(narray));
flightpathanglearray = totalmassarray;
groundtrackdistancearray = totalmassarray;

for i = 1:length(narray);
    n = narray(i);
    %Initial parameters
    m0 = 9400;
    h0 = 15000;    
    M0 = 5;
    rho0 = aer.getCheapDensity(h0);
    V0 = aer.getCheapSpeedofsound(h0)*M0;
    q0 = 0.5*rho0*V0^2;

    %Energy calculation
    Epot0 = m0*aer.getg(0)*h0;
    Ekin0 = 0.5*m0*V0^2;
    Etot0 = Epot0+Ekin0;

    %Velocity when no parachutes are used
    Vnodeceleration = sqrt(2*Etot0/m0);

    %Just retropropulsion
    Isp = 451;
    g0 = 9.81;
    Ve = Isp * g0;
    fuelfraction = 1-exp(-Vnodeceleration/Ve);
    fuelmass = m0*fuelfraction;

    R = 6;
    CD = 1.307;
    CDA = CD*pi*R^2;
    D0 = 0.5*rho0*V0^2*CDA;
    a0 = D0/m0;
    a_g0 = a0/g0;

    deceleration_gload = n;
    deceleration = deceleration_gload * g0;
    t_descent = V0/deceleration;
    d_3g = V0*t_descent-0.5*deceleration*t_descent^2;
    flightpathangle = acosd(h0/d_3g);
    groundtrackdistance = tand(flightpathangle)*h0;

    dt = 0.1;
    t = 0:dt:t_descent;
    h = h0*(1-t/t_descent);
    rho = aer.getCheapDensity(h);
    V = V0*(1-t/t_descent);
    q = 0.5*rho.*V.^2;
    CDA_parachute = 0.0*pi*15^2;
    D = q*CDA;
    D = D+(q*CDA_parachute).*(((D+q*CDA_parachute)/m0) < deceleration);
    a_gnothrust = D/m0/g0;
    F_thrust = (deceleration*m0-D+m0*aer.getg(h));
    F_thrust = F_thrust .* (F_thrust>0);
    mdot = F_thrust/Ve;
    
    M = V./aer.getCheapSpeedofsound(h);
    
    figure;
    hold on;
    skipframes = 30;
    plot(t(1:skipframes:end), (n*g0*m0+m0*aer.getg(h(1:skipframes:end)))/1000, '-o');
    plot(t(1:skipframes:end), D(1:skipframes:end)/1000, '-d');
    plot(t(1:skipframes:end), F_thrust(1:skipframes:end)/1000, '-s');
    xlabel('$t [s]$', 'interpreter', 'latex');
    ylabel('$F [kN]$', 'interpreter', 'latex');
    legend('Required force', 'Drag', 'Thrust', 'Location', 'east');
    axis([0, max(t), 0, 1.05*(max([D, F_thrust])/1000)]);
    grid on;
    matlab2tikz('plots/TDthrust.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
    
%     figure;
%     plot(t, M);
%     ylabel('M');
%     figure;
%     plot(t,a_gnothrust);
%     ylabel('g-force without thrust');

    thrustermass = sum(dt*mdot);
    thrusterdensity = 1; %kg/m3
    thrustervolume = thrustermass/thrusterdensity; %m^3
    tankmass = 2.7086e-8*thrustervolume^3-6.1703e-5*thrustervolume^2+6.629e-2*thrustervolume+1.3192;
    enginemass = 0.00144*max(F_thrust)+49.6;
    
    totalmassarray(i) = thrustermass+enginemass+tankmass;
    totalmassinteraction = thrustermass/2+enginemass+tankmass;
    groundtrackdistancearray(i) = groundtrackdistance;
    flightpathanglearray(i) = flightpathangle;
    
end



% figure;
% plot(narray, totalmassarray);
% xlabel('g_load');
% ylabel('mp');
% hold on;
% figure;
% plot(narray, flightpathanglearray);
% xlabel('g_load');
% ylabel('flight path angle');
% 
% figure;
% plot(narray, groundtrackdistancearray);
% xlabel('g_load');
% ylabel('d ground');

disp(' ');
disp(' ');
disp(strcat('Mach at start (10km): ', num2str(M0)));
disp(strcat('Velocity at start (10km): ', num2str(V0)));
disp(strcat('Velocity at surface without deceleration:', num2str(Vnodeceleration)));
disp(strcat('Fuel mass fraction when only retropropulsion is used:', num2str(fuelfraction)));
disp(strcat('Total fuel mass when only retropropulsion is used:', num2str(fuelmass)));
disp(strcat('Dynamic pressure at start (10km):', num2str(q0)));
disp(strcat('Flight path angle required:', num2str(flightpathangle)));
disp(strcat('Fuel mass (kg):', num2str(thrustermass)));
disp(strcat('Fuel mass/2 (interaction) (kg):', num2str(thrustermass/2)));
disp(strcat('Engine mass (kg):', num2str(enginemass)));
disp(strcat('Tank mass (kg):', num2str(tankmass)));
disp(strcat('Total mass, excluding parachute (kg):', num2str(totalmassarray(end))));
disp(strcat('Total mass/2 (interaction) (kg):', num2str(totalmassinteraction)));
disp(strcat('Deceleration at start (g):', num2str(a_g0)));
disp(strcat('Ground track distance (km):', num2str(groundtrackdistance/1000)));
disp(strcat('Maximum thrust delivered by the engine (N):', num2str(max(F_thrust))));



% h = 0:10000;
% speedofsound = aer.getCheapSpeedofsound(h);
% plot(h, speedofsound);
% % matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
