addpath('../matlab2tikz');

aer = marsatmosphere();


%Initial parameters
m0 = 10000;
h0 = 10000;
M0 = 4.55;
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
Ve = 4400;
fuelfraction = 1-exp(-Vnodeceleration/Ve);
fuelmass = m0*fuelfraction;

disp(strcat('Mach at start (10km): ', num2str(M0)));
disp(strcat('Velocity at start (10km): ', num2str(V0)));
disp(strcat('Velocity at surface without deceleration:', num2str(Vnodeceleration)));
disp(strcat('Fuel mass fraction when only retropropulsion is used:', num2str(fuelfraction)));
disp(strcat('Total fuel mass when only retropropulsion is used:', num2str(fuelmass)));
disp(strcat('Dynamic pressure at start (10km):', num2str(q0)));

R = 2.5;
CDA = pi*R^2;
D = 0.5*rho0*V0^2*CDA;
a = D/m0;
a_g = a/9.81;

deceleration_gload = 5;
deceleration = deceleration_gload * 9.81;
t_descent = V0/deceleration;
d_3g = V0*t_descent-0.5*deceleration*t_descent^2;
flightpathangle = acosd(h0/d_3g)

% h = 0:10000;
% speedofsound = aer.getCheapSpeedofsound(h);
% plot(h, speedofsound);
% % matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
