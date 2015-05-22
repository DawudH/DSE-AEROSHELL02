dt = [0.01,0.1,0.2,0.5,1];
e = [5.5, 6.5, 7.3, 10.0, 13.7];
cc = parula(3);
figure('name','Error between Kepler and numerical simulation');
h = loglog(dt,e,'d-');
grid on
set(h,'color',cc(1,:))
ylim([0 20])
xlabel('$\Delta t$ $\left[s\right]$','interpreter','latex','fontsize',15)
ylabel('Error in $R$ $\left[m\right]$','interpreter','latex','fontsize',15)
matlab2tikz('.\LaTeX\kep_num.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
