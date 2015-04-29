clc
clear all 
close all


m = 10000; %[kg]
CD = 1;
d = 12; % [m] diameter heatshield
S = 12^2*pi/4;
R_m = 6794000/2; %[m]
ry = 10* R_m; %[m]
v = 7000; %[m/s]
dt = 1;
h_atmos = 104 *10^3; % [m]
M_mars = 6.419*10^23; %[kg]
G = 6.673*10^-11; %[N*(m/kg)^2]
rx = -1.22*R_m:-0.001*R_m:-1.23*R_m;

cc = parula(length(rx)+3);
figure('name','orbits')
for i=1:length(rx)
    [out, R, V, A] = orbitmodel_new(rx(i),ry,R_m,m,CD,S,v,dt,h_atmos,M_mars,G);
    out
    t = 0:dt:(length(R)*dt-dt);
    subplot(3,1,1)
    hold on
    grid on
    plot(t,sqrt(R(:,1).^2 + R(:,2).^2 + R(:,3).^2),'color',cc(i,:));
    subplot(3,1,2)
    hold on
    grid on
    plot(t,sqrt(V(:,1).^2 + V(:,2).^2 + V(:,3).^2),'color',cc(i,:));
    subplot(3,1,3)
    hold on
    grid on
    plot(t,sqrt(A(:,1).^2 + A(:,2).^2 + A(:,3).^2),'color',cc(i,:));
end













% % circle plot:
% theta_plot = 0:0.01:2*pi;
% radius_mars = ones(1,length(theta_plot)) * R_m;
% radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + 150000);
% figure('name','Orbit')
% grid on
% axis equal
% hold on
% plot(R(:,1),R(:,2))
% polar(theta_plot,radius_mars,'r');
% polar(theta_plot,radius_mars_atmos,'g')
% 
% 
% figure('name','parameters over time')
% subplot(3,1,1)
% Rm = sqrt(R(:,1).^2 + R(:,3).^2 + R(:,2).^2);
% plot(t,Rm)
% grid on
% subplot(3,1,2)
% Vm = sqrt(V(:,1).^2 + V(:,3).^2 + V(:,2).^2);
% plot(t,Vm)
% grid on
% subplot(3,1,3)
% am = sqrt(a(:,1).^2 + a(:,3).^2 + a(:,2).^2);
% plot(t,am)
% grid on