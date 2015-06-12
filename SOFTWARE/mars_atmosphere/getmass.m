function y = getmass(x)
% close all;
% clear; clc;
aer = marsatmosphere();

g0 = x(1);
h0 = x(2);
M0 = x(3);

% for i = 1:length(garray);
    g = g0;
    %Initial parameters
    m0 = 10000;
%     h0 = 10000;
%     M0 = 3;
    V0 = aer.getCheapSpeedofsound(h0)*M0;


    R = 6;
    CD = 1.22;
    CDA = CD*pi*R^2;

    deceleration_gload = g;
    deceleration = deceleration_gload * 9.81;
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
    CDA_parachute = 0.3*pi*15^2;
    D = q*CDA;
    D = D+(q*CDA_parachute).*(((D+q*CDA_parachute)/m0) < deceleration);
    F_thrust = (deceleration*m0-D+m0*aer.getg(h));
    F_thrust = F_thrust .* (F_thrust>0);
    

    SFC = 0.225e-3; %kg/N/s http://en.wikipedia.org/wiki/Specific_impulse
    thrustermass = SFC*0.1*sum(F_thrust(F_thrust>0));
    enginemass = 0.00144*max(F_thrust)+49.6;
    
    y = [thrustermass+enginemass, flightpathangle];

    
end


% figure;
% plot(garray, thrustermassarray);
% xlabel('g_load');
% ylabel('mp');
% hold on;
% figure;
% plot(garray, flightpathanglearray);
% xlabel('g_load');
% ylabel('flight path angle');
% 
% figure;
% plot(garray, groundtrackdistancearray);
% xlabel('g_load');
% ylabel('d ground');

% disp(strcat('Mach at start (10km): ', num2str(M0)));
% disp(strcat('Velocity at start (10km): ', num2str(V0)));
% disp(strcat('Velocity at surface without deceleration:', num2str(Vnodeceleration)));
% disp(strcat('Fuel mass fraction when only retropropulsion is used:', num2str(fuelfraction)));
% disp(strcat('Total fuel mass when only retropropulsion is used:', num2str(fuelmass)));
% disp(strcat('Dynamic pressure at start (10km):', num2str(q0)));
% disp(strcat('Flight path angle required:', num2str(flightpathangle)));
% disp(strcat('Fuel mass (kg):', num2str(thrustermass)));
% disp(strcat('Engine mass (kg):', num2str(enginemass)));
% disp(strcat('Total mass, excluding parachute (kg):', num2str(totalmassarray(end))));
% disp(strcat('Deceleration at start (g):', num2str(a_g0)));
% disp(strcat('Ground track distance (km):', num2str(groundtrackdistance/1000)));



% h = 0:10000;
% speedofsound = aer.getCheapSpeedofsound(h);
% plot(h, speedofsound);
% % matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
