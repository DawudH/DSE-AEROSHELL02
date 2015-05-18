clear
close all
clc

%% Do the check

out = discretization_check([0.00625 0.0125 0.025 0.05 0.1 0.2 0.4 0.8]);

%% Plot the results

figure('name','discretization verification')
cc = parula(3);
h = loglog(out.dt,out.max_error,'o-');
set(h,'color',cc(2,:));
grid on
xlabel('$\Delta t$ $\left[s\right]$','interpreter','latex','fontsize',15)
ylabel('max error in $R$ $\left[m\right]$','interpreter','latex','fontsize',15)