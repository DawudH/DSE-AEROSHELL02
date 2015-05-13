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
%R wrt the elipse reference frame
R = r*[cos(theta),sin(theta),0]
%Express in 0-reference frame
R0 = rotz(rad2deg(theta_p))*R';
x = R0(1);

%Velocity
%magnitude
v = sqrt(mu*(2/r-1/a));
%direction
x_rc = (x+a*e);
rc = b^2*x_rc/(a^2*sqrt(b^2*(x^2-a^2)/a^2));
V_unit = [1,rc,0]/norm([1,rc,0])
V = v*V_unit;

%Acceleration
A = -mu/r^3*R

%%Output
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