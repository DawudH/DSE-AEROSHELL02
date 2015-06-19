clear
close all
clc

addpath('..\aerodynamic_coefficients')
addpath('..\matlab2tikz')
aero_coef = aeroProperties();
alpha = [0:30]*pi / 180;

CXA = aero_coef.getCXA(alpha);
CZA = aero_coef.getCZA(alpha);
[CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);


x = CXA;
y = 0;
z = CZA;

Mx = 0; % 1/1
My = CMYA; % 1/300
Mz = 0; % 1/3

dx = [-1,-3,-5];
cc = parula(4);
marker = ['-o','-d','-x'];
figure('name','dz over alpha')
hold on
grid on
for i = 1:length(dx)
    dz = (My+z.*dx(i))./x;
    dy = (Mz-y.*dx(i))./-x;
    plot(alpha*180/pi,dz,marker(2*i-1:2*i),'color',cc(i,:))
    legend_str{i} = ['Xcg = ' num2str(dx(i)) ' [m]'];
end
h=legend(legend_str,'Location','SouthWest');

xlabel('$\alpha$ $\left[deg\right]$','interpreter','latex','fontsize',15)
ylabel('$Z_{c.g.}$ $\left[m\right]$','interpreter','latex','fontsize',15)
matlab2tikz('LaTeX\moment.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);