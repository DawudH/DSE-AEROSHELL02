addpath('../matlab2tikz');
clear; close all

x = 0:0.1:pi/2;
L = cos(x).^2.*sin(x);
D = cos(x).^3;
LD = L./D;

figure
grid on
hold on
axis([0 rad2deg(pi/2) 0 2])
ylabel('Relative value', 'interpreter', 'latex')
xlabel('Incidence angle [deg]', 'interpreter', 'latex')

plot(rad2deg(x),L,'-o')
plot(rad2deg(x),D,'-x')
plot(rad2deg(x),LD,'-+')

h = legend('Relative Lift', 'Relative Drag', 'Relative Lift/Drag');
h.Interpreter = 'latex';
matlab2tikz('plots/LDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

