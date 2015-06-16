close all;clear
addpath('../matlab2tikz')
figure
grid on

cc = parula(5);
x = [0 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
hold on
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-o','color',cc(1,:));

x = [0.5 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-x','color',cc(2,:));

x = [1 2 1 0 0];
[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
plot(rad2deg(mod.alpha_array), mod.CLCD_array, '-+','color',cc(3,:));

xlabel('Angle of attack ($\alpha$) [deg]', 'interpreter', 'latex');
ylabel('$C_{L}C^{-1}_{D}$','interpreter', 'latex');
h = legend('No offset', '0.5m offset' , '1m offset');
h.Interpreter = 'latex';
h.Location = 'SouthEast'


matlab2tikz('plots/LoverDplot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);