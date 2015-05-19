clear
close all;

aerirve = aeroProperties('irve');
aerapollo = aeroProperties('apollo');
aerisotensoid = aeroProperties('isotensoid');
aertorus = aeroProperties('torus');

alpha = 0:0.01:deg2rad(60);

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCLA(alpha));
plot(rad2deg(alpha), aerapollo.getCLA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCLA(alpha));
plot(rad2deg(alpha), aertorus.getCLA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_LA [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/cl', '-depsc');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCDA(alpha));
plot(rad2deg(alpha), aerapollo.getCDA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCDA(alpha));
plot(rad2deg(alpha), aertorus.getCDA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_DA [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/cd', '-depsc');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCMYA(alpha));
plot(rad2deg(alpha), aerapollo.getCMYA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCMYA(alpha));
plot(rad2deg(alpha), aertorus.getCMYA(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{M}A [m^2]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/cm', '-depsc');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getLiftGradient(alpha));
plot(rad2deg(alpha), aerapollo.getLiftGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getLiftGradient(alpha));
plot(rad2deg(alpha), aertorus.getLiftGradient(alpha));
xlabel('$\alpha [deg]$');
ylabel('$C_{L_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/clalpha', '-depsc');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getDragGradient(alpha));
plot(rad2deg(alpha), aerapollo.getDragGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getDragGradient(alpha));
plot(rad2deg(alpha), aertorus.getDragGradient(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{D_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/cdalpha', '-depsc');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getMomentGradient(alpha));
plot(rad2deg(alpha), aerapollo.getMomentGradient(alpha));
plot(rad2deg(alpha), aerisotensoid.getMomentGradient(alpha));
plot(rad2deg(alpha), aertorus.getMomentGradient(alpha));
xlabel('$\alpha [deg]$', 'interpreter', 'latex');
ylabel('$C_{M_\alpha}A [\frac{m^2}{rad}]$', 'interpreter', 'latex');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
print('plots/cmalpha', '-depsc');
