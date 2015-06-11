addpath('../matlab2tikz');

aer = marsatmosphere();


%Initial parameters
m0 = 10000;
h0 = 10000;
M0 = 5;
V0 = aer.getCheapSpeedofsound(0)*M0;
q0 = 0.5*aer.getCheapDensity(h0)*V0^2;

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

% h = 0:10000;
% speedofsound = aer.getCheapSpeedofsound(h);
% plot(h, speedofsound);
% % matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
