load('effect_mu.mat')
results.gamma = 21.83:0.001:21.865;
figure('name','param vs. mu')
subplot(1,3,1)
hold on
grid on
ylabel('$q$ $\left[Pa\right]$','interpreter','latex')
plot(results.mu,results.q)
subplot(1,3,2)
hold on
grid on
ylabel('$h$ $\left[m\right]$','interpreter','latex')
plot(results.mu,results.h)
subplot(1,3,3)
hold on
grid on
%xlim([21.83,21.87])
% ylim([41.68,41.69])
ylabel('$M$ $\left[-\right]$','interpreter','latex')
plot(results.mu,results.M)