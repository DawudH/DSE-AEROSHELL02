clear
close all
clc

addpath('..\aerodynamic_coefficients')
aero_coef = aeroProperties('irve');
alpha = [0:20]*pi / 180;

[CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);

x = CDA;
y = 0;
z = CLA;

Mx = 0; % 1/1
My = CMYA; % 1/300
Mz = 0; % 1/3

dx = [1,3,5];
figure('name','dz over alpha')
hold on
for i = 1:length(dx)

    dz = (My+z.*dx(i))./x;
    dy = (Mz-y.*dx(i))./-x;
    plot(alpha*180/pi,dz)
    legend_str{i} = ['Xcg = ' num2str(dx(i)) ' [m]'];
end
grid on
h=legend(legend_str);

xlabel('$\alpha$ $\left[^\circ\right]$','interpreter','latex','fontsize',15)
ylabel('$Z_{c.g.}$ $\left[m\right]$','interpreter','latex','fontsize',15)
