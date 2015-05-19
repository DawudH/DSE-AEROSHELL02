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
xlabel('\alpha (deg)');
ylabel('C_LA (m^2)');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCDA(alpha));
plot(rad2deg(alpha), aerapollo.getCDA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCDA(alpha));
plot(rad2deg(alpha), aertorus.getCDA(alpha));
xlabel('\alpha (deg)');
ylabel('C_DA (m^2)');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');

figure;
hold on;
plot(rad2deg(alpha), aerirve.getCMYA(alpha));
plot(rad2deg(alpha), aerapollo.getCMYA(alpha));
plot(rad2deg(alpha), aerisotensoid.getCMYA(alpha));
plot(rad2deg(alpha), aertorus.getCMYA(alpha));
xlabel('\alpha (deg)');
ylabel('C_{M_Y}A (m^2)');
legend('Stacked Toroid, Tension Cone', 'Rigid', 'Isotensoid', 'Trailing Ballute');
