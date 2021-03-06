clear
close all
clc

added_paths

au = 1.496e11; %m
r_e = 1.013 * au; %m
r_m = 1.53 * au; %m
R_m_SOI = 5.77e8;
R_earth = 6371e3; % radius surface of the earth [m]
M_earth = 5.97219e24; % mass earth
M_mars = 0.64174*10^24; %[kg]
G = 6.67384*10^-11;
M_sun = 1.988435e30; % kg
R_m = 3.389945945211271e6; %[m]
h_atm = 400e3; %[m] Height of atmosphere

V_e = sqrt(G * M_sun / r_e);
V_m = sqrt(G * M_sun / r_m);
V_esc_earth = sqrt(2*G*M_earth / R_earth);

delta_V_hohmann = sqrt( G*M_sun * (2/r_e - 2 / (r_m + r_e))) - V_e;
delta_V = delta_V_hohmann:5e1:delta_V_hohmann+6300;

V_sc_p = delta_V + V_e;

a =  (2 / r_e - V_sc_p.^2 / ( G * M_sun)).^(-1);
e = 1- r_e ./ a;
b = sqrt(a.^2 - (e.*a).^2);

V_sc_m = sqrt( G * M_sun * (2/r_m - 1./a));
a_m = (2 / R_m_SOI - (V_sc_m-V_m).^2 / ( G * M_mars)).^(-1);
V_sc_m_atmos = sqrt(G*M_mars*(2/(R_m + h_atm) - 1./a_m));

theta = acos( (a .* (1-e.^2) / r_m - 1) ./ e);
A = zeros(1,length(delta_V));
for i = 1:length(delta_V)
    fun = @(TH) 1/2 * (a(i) * (1-e(i)^2) ./ (1 + e(i)*cos(TH))).^2;
    A(i) = integral(fun,0,theta(i));
end

t = 2*A ./ sqrt(a*G*M_sun.*(1-e.^2))/3600/24;
cc = parula(6);
index = V_sc_m - V_m >= 0;
V_sc_m_atmos_e = V_sc_m_atmos(index);
delta_V_e = delta_V(index);
figure('name','delta_v versus transfer time and SOI of Mars')
hold on
marker_index = 1:5:length(delta_V);
marker_index2 = 1:5:length(delta_V_e);
[AX,H1,H2] = plotyy((delta_V(marker_index)+V_esc_earth)/1000,t(marker_index),(delta_V_e(marker_index2)+V_esc_earth)/1000,V_sc_m_atmos_e(marker_index2)/1000);
H1.Color = cc(1,:);
H1.Marker = 'o';
H2.Color = cc(3,:);
H2.Marker = '*';
H1.LineStyle = 'none';
H2.LineStyle = 'none';
AX(1).YColor = cc(1,:);
AX(2).YColor = cc(3,:);
[AX,H1,H2] = plotyy((delta_V+V_esc_earth)/1000,t,(delta_V(index)+V_esc_earth)/1000,V_sc_m_atmos(index)/1000);
H1.Color = cc(1,:);
H2.Color = cc(3,:);
AX(1).YColor = cc(1,:);
AX(2).YColor = cc(3,:);
AX(1).YLabel.Interpreter = 'latex';
AX(1).XLabel.Interpreter = 'latex';
AX(2).YLabel.Interpreter = 'latex';
AX(2).YLabel.String = '$V_{s/c}$ $\left[ km \cdot s^{-1}\right]$';
xlabel('$\Delta V$ from earth $\left[ km \cdot s^{-1} \right]$')
ylabel('Transfer time $\left[ days \right]$')
h5 = line([1.9623e+04 1.9623e+04]/1000, ylim,'linestyle','-.','color',cc(4,:));
legend([H1 H2 h5],'Transfer time',...
        'Vs/c at boundary of atmosphere of Mars',...
        '7 km/s at atmosphere boundary',...
        'location','northWest')
grid on
%matlab2tikz('.\LaTeX\transfer_time.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure('name','transfer orbit')
hold on
axis equal
theta_plot = 0:0.01:2*pi;
index_orbit = 112;
theta_mars = acos(a(index_orbit)* (1-e(index_orbit)^2) / (e(index_orbit) * r_m) - 1/e(index_orbit));
theta_plot_sc = 0:0.01:theta_mars;
r_sc = a(index_orbit) * (1-e(index_orbit)^2) ./ (1 + e(index_orbit) * cos(theta_plot_sc) );
r_mars = ones(1,length(theta_plot)) * r_m;
r_earth = ones(1,length(theta_plot)) * r_e;
R_sun = ones(1,length(theta_plot)) * 10e9;
plot(0,0,'o','markerfacecolor',[1 0.8 0],'markersize',19,'markeredgecolor','none')
plot(r_e,0,'o','markerfacecolor',[30 144 255]/255,'markersize',11,'markeredgecolor','none')
plot(cos(theta_mars)*r_m,sin(theta_mars)*r_m,'o','markerfacecolor',[161 37 27]/255,'markersize',9,'markeredgecolor','none')
hre = polar(theta_plot,r_earth,'--'); 
hre.Color = [30 144 255]/255;
hrm = polar(theta_plot,r_mars,'--');
hrm.Color = [161 37 27]/255;
hsc = polar(theta_plot_sc,r_sc,'-.'); 
hsc.Color = 'k';
set(gca,'Visible','off')
set(gcf,'color',[1 1 1])
legend('Sun','Earth','Mars','Earth trajectory','Mars trajectory','Interplanetary trajectory','location','westoutside')
% matlab2tikz('.\LaTeX\transfer_orbit.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

figure('name','transfer orbit black back')
hold on
axis equal
theta_plot = 0:0.01:2*pi;
index_orbit = 112;
theta_mars = acos(a(index_orbit)* (1-e(index_orbit)^2) / (e(index_orbit) * r_m) - 1/e(index_orbit));
theta_plot_sc = 0:0.01:theta_mars;
r_sc = a(index_orbit) * (1-e(index_orbit)^2) ./ (1 + e(index_orbit) * cos(theta_plot_sc) );
r_mars = ones(1,length(theta_plot)) * r_m;
r_earth = ones(1,length(theta_plot)) * r_e;
R_sun = ones(1,length(theta_plot)) * 10e9;
plot(0,0,'o','markerfacecolor',[1 0.8 0],'markersize',19,'markeredgecolor','none')
plot(r_e,0,'o','markerfacecolor',[30 144 255]/255,'markersize',11,'markeredgecolor','none')
plot(cos(theta_mars)*r_m,sin(theta_mars)*r_m,'o','markerfacecolor',[229 52 39]/255,'markersize',9,'markeredgecolor','none')
hsc = polar(theta_plot_sc,r_sc,'-.'); 
hsc.Color = 'w';
hre = polar(theta_plot,r_earth,'--'); 
hre.Color = [30 144 255]/255;
hrm = polar(theta_plot,r_mars,'--');
hrm.Color = [229 52 39]/255;
set(gca,'Visible','off')
set(gcf,'color',[0 0 0])
matlab2tikz('.\LaTeX\transfer_orbit_black.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);

