close all;clear
addpath('../matlab2tikz')
figure
grid on

x = [0 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
hold on
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-o');

x = [0.5 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-x');

x = [1 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-+');

xlabel('Angle of attack ($\alpha$) [deg]', 'interpreter', 'latex');
ylabel('$\frac{C_{L}}{C_{D}}$','interpreter', 'latex');
h = legend('No offset', '0.5m offset' , '1m offset');
h.Interpreter = 'latex';
h.Location = 'SouthEast'


matlab2tikz('plots/LoverDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);