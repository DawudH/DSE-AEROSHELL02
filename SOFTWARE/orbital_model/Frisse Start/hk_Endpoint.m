function [ out ] = hk_Endpoint(param,G,M,R_m,h_atm)
%Calculates the endpoint of the hyperbolic entry orbit
%%input
a = param.a;
e = param.e;
b = param.b;
theta_p = param.theta_p;
mu = G*M;

%%Calculation

%Position
r  = R_m+h_atm;
theta = -acos((a*(1-e^2)-r)/(r*e));
theta0 = theta_p -2*pi + theta;
%R wrt the elipse reference frame
R0 = r*[cos(theta0),sin(theta0),0];
x = R0(1);
%Velocity
%magnitude
v = sqrt(mu*(2/r-1/a));
%direction
rc = b^2*(a*e+x)/(a^2*sqrt(b^2*(a*e+x)^2/a^2-b^2));
%rc = b^2*x_rc/(a^2*sqrt(b^2*(x_rc^2-a^2)/a^2))
%rc = a^2*x_rc/(b^2*sqrt(a^2*(b^2+x_rc^2)/b^2))
V_unit = [1,rc,0]/norm([1,rc,0]);
V_unit = (rotz(-theta_p)*V_unit')';
V = v*V_unit;

%Acceleration
A = -mu/r^3*R0;

%%Output
out.rc = rc;
out.R = R0;
out.speed_sound = 0;
out.V = V;
out.M = 0;
out.A = A;
out.Ag = A;
out.Ad = [0,0,0];
out.Al = [0,0,0];
out.J = [0,0,0];
out.q = 0;
out.theta = theta;
end