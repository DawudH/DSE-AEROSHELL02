load('effect_mu.mat')
results.gamma = 21.83:0.001:21.865;
figure('name','param vs. mu')
subplot(1,3,1)
hold on
grid on
ylabel('$q$ $\left[Pa\right]$','interpreter','latex')
xlabel('$\mu$ $\left[deg\right]$','interpreter','latex')
plot(results.mu,results.q)
subplot(1,3,2)
hold on
grid on
xlabel('$\mu$ $\left[deg\right]$','interpreter','latex')
ylabel('$h$ $\left[km\right]$','interpreter','latex')
plot(results.mu,results.h/1000)
subplot(1,3,3)
hold on
grid on
ylim([41.68,41.69])
xlabel('$\mu$ $\left[deg\right]$','interpreter','latex')
ylabel('$M$ $\left[-\right]$','interpreter','latex')
plot(results.mu,results.M)

matlab2tikz('.\LaTeX\effectmu.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
