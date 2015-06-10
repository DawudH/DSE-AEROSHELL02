addpath('../matlab2tikz');
clear; close all

x = 0:0.05:pi/2;
L = cos(x).^2.*sin(x);
D = cos(x).^3;
LD = L./D;

figure
hold on
axis([0 rad2deg(pi/2) 0 2])
ylabel('Relative value')
xlabel('Incidence angle [deg]', 'interpreter', 'latex')

plot(rad2deg(x),L)
plot(rad2deg(x),D)
plot(rad2deg(x),LD)

legend('Relative Lift', 'Relative Drag', 'Relative Lift/Drag')

matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

