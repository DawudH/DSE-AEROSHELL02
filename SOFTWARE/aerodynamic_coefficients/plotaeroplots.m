clear
close all;
addpath('..\matlab2tikz');

aerirve = aeroProperties('irve');
aerapollo = aeroProperties('apollo');
aerisotensoid = aeroProperties('isotensoid');
aerballute = aeroProperties('ballute');

alpha = linspace(0, deg2rad(60), 60);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCLA(alpha));
plot(rad2deg(alpha), aerapollo.getCLA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCLA(alpha));
plot(rad2deg(alpha), aerballute.getCLA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_LA [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cl.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCDA(alpha));
plot(rad2deg(alpha), aerapollo.getCDA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCDA(alpha));
plot(rad2deg(alpha), aerballute.getCDA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_DA [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cd.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCMYA(alpha));
plot(rad2deg(alpha), aerapollo.getCMYA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCMYA(alpha));
plot(rad2deg(alpha), aerballute.getCMYA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{M}A [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cm.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getLiftGradient(alpha));
plot(rad2deg(alpha), aerapollo.getLiftGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getLiftGradient(alpha));
plot(rad2deg(alpha), aerballute.getLiftGradient(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{L_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\clalpha.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getDragGradient(alpha));
plot(rad2deg(alpha), aerapollo.getDragGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getDragGradient(alpha));
plot(rad2deg(alpha), aerballute.getDragGradient(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{D_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cdalpha.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getMomentGradient(alpha));
plot(rad2deg(alpha), aerapollo.getMomentGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getMomentGradient(alpha));
plot(rad2deg(alpha), aerballute.getMomentGradient(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{M_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cmalpha.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCLCD(alpha));
plot(rad2deg(alpha), aerapollo.getCLCD(alpha));
plot(rad2deg(alpha), aerisotensoid.getCLCD(alpha));
plot(rad2deg(alpha), aerballute.getCLCD(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$\frac{C_L}{C_D} [-]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\clcd.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);


figure;
hold on;
alpha = alpha(2:end);
plot(rad2deg(alpha), aerirve.getCMCL(alpha));
plot(rad2deg(alpha), aerapollo.getCMCL(alpha));
plot(rad2deg(alpha), aerisotensoid.getCMCL(alpha));
plot(rad2deg(alpha), aerballute.getCMCL(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$\frac{C_M}{C_L} [-]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute', 'Location', 'northoutside');
grid on;
matlab2tikz('.\plots\cmcl.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

