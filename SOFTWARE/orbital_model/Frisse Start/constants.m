%%Universal constants
%Gravitational constant
G = 6.67384*10^-11; %[N*(m/kg)^2]

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
M_mars = 0.64174*10^24; %[kg]
%Radius of Mars

R_m = 3.389945945211271e6; %[m]
%Height of atmosphere
h_atm = 400e3; %[m]
%Sphere of influence
SOI = 0.576e9; %[m]

% crash margin above the surface
crash_margin = 15000; %[m]


%Earth gravitational accaleration
g_earth = 9.81; %[m/s^2]

%Height of end of mission
height_EOM = 15000; %[m]

