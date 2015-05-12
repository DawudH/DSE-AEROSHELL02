addpath('..\matlab2tikz')
clc
clear
close all

figure(1)
Z = peaks(20);
contourf(Z,10)
grid on
xl = xlim;
yl = ylim;
matlab2tikz('.\LaTeX\test.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
fix_grid_contourf('.\LaTeX\test.tikz')

figure(2)
plot([min(xl)-1,min(xl)-1],[min(yl)-1,min(yl)-1]);
grid on
ylim(yl)
xlim(xl)
matlab2tikz('.\LaTeX\grid.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);