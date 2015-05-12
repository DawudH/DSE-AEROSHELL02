function [out] = hk_Parameters(R0,V0,G,M)
%Calculates the kepler orbit parameters needed to calculate the hyperbolic
%flightpath

%%Input
r0 = norm(R0);
v0 = norm(V0);
theta0 = atan2(R0(2),R0(1));
mu = G*M;

%%Calculations
a = (2/r0-v0^2/mu)^-1;
e = -r0*cos(theta0)/(2*a)+sqrt(r0^2*cos(theta0)^2/(4*a^2)-r0/a+1);
b = a*sqrt(e^2-1);
%%Output
out.a = a;
out.e = e;
out.b = b;
out
end

