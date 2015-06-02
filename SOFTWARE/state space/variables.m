%%%constant
%%Sebastiaan
g0 = 3.75; % [m/s^2]

%%Twan
m = 10000; %[kg]
Sref = pi/4*12^2; %[m^2]


%%%variable
%%Twan
V0 = 7000; %[m/s]
M0 = 50; %[-]
rho0 = 1e-15; %[kg/m^3]
q0 = rho0*V0^2/2; %[kg/(m*s^2)]
CD0 = ??; %[-]
D0 = CD0*q0*Sref;

%%Sebastiaan
R0 = 400*10^3 + 3.389945945211271e6; % [m]
CL0 = getCL(alpha0); % [-]
L0 = CL0*q0*Sref;
sigma0 = 0; % [-]