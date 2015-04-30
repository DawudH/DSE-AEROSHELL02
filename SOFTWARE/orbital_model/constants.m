 omega_m = 2*pi / (24.6229622 * 3600 ); %[rad/sec]
Omega_m = [0,0,1]*omega_m; %[rad/sec]
m = 10000; %[kg]
d = 12; % [m] diameter heatshield
S = 12^2*pi/4; %[m^2]
R_m = 6794000/2; %[m]
h_atm = 104e3; %[m]
z = [0 0 1]; %[-]

%save('constants.mat');