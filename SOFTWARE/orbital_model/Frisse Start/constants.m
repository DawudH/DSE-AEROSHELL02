addpath('..\..\mars_atmosphere')

%%Universal constants
%Gravitational constant
G = 6.673*10^-11; %[N*(m/kg)^2]

%%Geometric properties of s/c
%Mass
m = 10000; %[kg]
%Diameter of the heatshield
d = 12; %[m]
%Area of the heatshield
S = d^2*pi/4; %[m^2]

%%Mars parameters
%Rotation of mars and its atmosphere
omega_m = 2*pi / (24.6229622 * 3600 ); %[rad/sec]
Omega_m = [0,0,1]*omega_m; %[rad/sec]
%Mass
M_mars = 6.419*10^23; %[kg]
%Radius of Mars
R_m = 6794000/2; %[m]
%Height of atmosphere
h_atm = 400e3; %[m]
%Sphere of influence
SOI = 0.576e9; %[m]


%Earth gravitational accaleration
g_earth = 9.81; %[m/s^2]

%Height of end of mission
height_EOM = 10000; %[m]

