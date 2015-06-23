clc
clear
close all

addpath('.\..\matlab2tikz')

% generate atmos model output plots

atm = marsatmosphere(); 

lat = 0;
lon = 0:1:350;
lat2 = -80:1:80;
lon_fixed = 180;

h = 0:1000:400e3;
h_fixed = 200e3;
h_fixed2 = 100e3;

rho = atm.getDensity(lat,lon,h_fixed);
T = atm.getTemperature(lat,lon,h_fixed);

fontsize = 15;

cc = parula(4);
figure('name', ['Atmospheric model dependency on longitude for a latitude of 0 deg, at ' num2str(h_fixed/1000) ' [km] high'])
hold on
[hAx,hrho,hT] = plotyy(lon,rho,lon,T);
xlim(hAx(1),[0 360])
xlim(hAx(2),[0 360])
set(hAx(1),'Ytick',linspace(min(ylim(hAx(1))),max(ylim(hAx(1))),5))
set(hAx(2),'Ytick',linspace(min(ylim(hAx(2))),max(ylim(hAx(2))),5))
set(hAx(2),'XTick',[]);
set(hrho,'color',cc(1,:))
set(hAx(1),'Ycolor',cc(1,:))
set(hT,'color',cc(2,:))
set(hAx(2),'Ycolor',cc(2,:))
xlim(hAx(1),[0 360])
xlim(hAx(2),[0 360])
grid on
xlabel('Longitude $\left[^\circ\right]$','interpreter','latex','fontsize',fontsize)
ylabel(hAx(1),'Density $\left[\frac{kg}{m^3}\right]$','interpreter','latex','fontsize',fontsize)
ylabel(hAx(2),'Temperature $\left[K\right]$','interpreter','latex','fontsize',fontsize)
matlab2tikz('.\LaTeX\lon_25.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

rho = atm.getDensity(lat,lon,h_fixed2);
T = atm.getTemperature(lat,lon,h_fixed2);

cc = parula(4);
figure('name', ['Atmospheric model dependency on longitude for a latitude of 0 deg, at ' num2str(h_fixed2/1000) ' [km] high'])
hold on
[hAx,hrho,hT] = plotyy(lon,rho,lon,T);
xlim(hAx(1),[0 360])
xlim(hAx(2),[0 360])
set(hAx(1),'Ytick',linspace(min(ylim(hAx(1))),max(ylim(hAx(1))),5))
set(hAx(2),'Ytick',linspace(min(ylim(hAx(2))),max(ylim(hAx(2))),5))
set(hAx(2),'XTick',[]);
set(hrho,'color',cc(1,:))
set(hAx(1),'Ycolor',cc(1,:))
set(hT,'color',cc(2,:))
set(hAx(2),'Ycolor',cc(2,:))
grid on
xlabel('Longitude $\left[^\circ\right]$','interpreter','latex','fontsize',fontsize)
ylabel(hAx(1),'Density $\left[\frac{kg}{m^3}\right]$','interpreter','latex','fontsize',fontsize)
ylabel(hAx(2),'Temperature $\left[K\right]$','interpreter','latex','fontsize',fontsize)
matlab2tikz('.\LaTeX\lon_50.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);



rho = zeros(length(h),1);
T = zeros(length(h),1);


for i = 1:length(h)
    rho(i) = atm.getDensity(lat,lon_fixed,h(i));
    T(i) = atm.getTemperature(lat,lon_fixed,h(i));
end

figure('name', 'Atmospheric density dependency on hight for a longitude of 180 a latitude of 0 deg')
semilogx(rho,h/1000,'color',cc(1,:));
ax = gca;
grid on
ylabel('Height $\left[km\right]$','interpreter','latex','fontsize',fontsize)
xlabel('Density $\left[\frac{kg}{m^3}\right]$','interpreter','latex','fontsize',fontsize)
matlab2tikz('.\LaTeX\rho_h_180_0.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure('name', 'Atmospheric temperature model dependency on hight for a longitude of 180 a latitude of 0 deg')
plot(T,h/1000,'color',cc(2,:))
xlim([115 210])
grid on
ylabel('Height $\left[km\right]$','interpreter','latex','fontsize',fontsize)
xlabel('Temperature $\left[K\right]$','interpreter','latex','fontsize',fontsize)
matlab2tikz('.\LaTeX\T_h_180_0.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

rho = atm.getDensity(lat2,lon,h_fixed);
[X,Y] = meshgrid(lon,lat2);
figure('name', 'Atmospheric density model dependency longitude and latitude')
contourf(X,Y,rho','ShowText','on')
grid on
ylabel('latitude $\left[^\circ\right]$','interpreter','latex','fontsize',fontsize)
xlabel('Longitude $\left[^\circ\right]$','interpreter','latex','fontsize',fontsize)
%matlab2tikz('.\LaTeX\T_h_180_0.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
