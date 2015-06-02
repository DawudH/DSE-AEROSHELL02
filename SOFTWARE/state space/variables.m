function [g0,m,bref,Sref,V0,M0,rho0,alpha0,gamma0,gammadot0,q0,CL0A,CD0A,CMY0A,D0,L0,MY0,R0,sigma0,Ixx,Iyy,Izz]=variables()
% create aerocoef object
aero_coef = aeroProperties('irve');

%%%constant
%%Sebastiaan
g0 = 3.75; % [m/s^2]

%%Twan
m = 10000; %[kg]
bref = 12; %[m]
Sref = (pi/4)*bref^2; %[m^2]


%%%variable
%%Twan
V0 = 7000; %[m/s]
M0 = 50; %[-]
rho0 = 1e-15; %[kg/m^3]
alpha0 = 0; %[rad]
%flight path angle
gamma0 = 5/180*pi; %[rad]
gammadot0 = 0/180*pi; %[rad/s]
%dynamic pressure
q0 = rho0*V0^2/2; %[kg/(m*s^2)]
%Aerodynamic forces
[CL0A, CD0A, CMY0A] = aero_coef.aeroCoeffs(alpha0);
D0 = CD0A*q0;
L0 = CL0A*q0;
MY0 = CMY0A*q0;

%%Sebastiaan
R0 = 400e3 + 3.389945945211271e6; % [m]
sigma0 = 0/180*pi; % [-]
Ixx = 10e3; % mass moment of inertia
Iyy = 10e3; % wild guesses
Izz = 10e3;
end