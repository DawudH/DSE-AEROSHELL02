load('effect_alpha.mat')
results.gamma = 21.83:0.001:21.865;
figure('name','param vs. alpha')
subplot(1,3,1)
hold on
grid on
%xlim([0 60])
%xlim([21.83 21.87])
ylabel('$q$ $\left[Pa\right]$','interpreter','latex')
xlabel('$\alpha$ $\left[deg\right]$','interpreter','latex')
plot(results.alpha,results.q)
subplot(1,3,2)
hold on
grid on
%xlim([0 60])
%xlim([21.83 21.87])
xlabel('$\alpha$ $\left[deg\right]$','interpreter','latex')
ylabel('$h$ $\left[km\right]$','interpreter','latex')
plot(results.alpha,results.h/1000)
subplot(1,3,3)
hold on
grid on
%xlim([0 60])
%xlim([21.83 21.87])
ylim([41.6,41.7])
xlabel('$\alpha$ $\left[deg\right]$','interpreter','latex')
ylabel('$M$ $\left[-\right]$','interpreter','latex')
plot(results.alpha,results.M)

matlab2tikz('.\LaTeX\effectalpha.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
