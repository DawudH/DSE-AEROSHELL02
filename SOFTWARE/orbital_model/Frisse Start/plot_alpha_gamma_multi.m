close all

load('phi_0.mat')
load('phi_30.mat')
load('phi_60.mat')


% extrapolated
figure('name','fitted to datapoints and extrapolated')
hold on
plot(21.65:0.001:22,polyval(phi_0.flyby.polyfit,21.65:0.001:22)*180/pi)
plot(21.85:0.001:22,polyval(phi_0.accel.polyfit,21.85:0.001:22)*180/pi)
plot(21.74:0.001:21.93,polyval(phi_0.crash.polyfit,21.74:0.001:21.93)*180/pi)

%plot(21.74:0.001:22,polyval(phi_30.flyby.polyfit,21.74:0.001:22)*180/pi)
%plot(21.85:0.001:22,polyval(phi_30.accel.polyfit,21.85:0.001:22)*180/pi)
%plot(21.74:0.001:21.93,polyval(phi_30.crash.polyfit,21.74:0.001:21.93)*180/pi)

%plot(21.74:0.001:22,polyval(phi_60.flyby.polyfit,21.74:0.001:22)*180/pi)
plot(21.83:0.001:22,polyval(phi_60.accel.polyfit,21.83:0.001:22)*180/pi)
plot(21.65:0.001:21.93,polyval(phi_60.crash.polyfit,21.65:0.001:21.93)*180/pi)
legend('phi 0: flyby','phi 0: 3g','phi 0: crash','phi 60: 3g','phi 60: crash')
ylim([-5 35])
grid on
