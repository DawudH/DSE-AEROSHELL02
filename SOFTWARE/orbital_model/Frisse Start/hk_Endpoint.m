function [ out ] = hk_Endpoint(param,G,M,R_m,h_atm)
%Calculates the endpoint of the hyperbolic entry orbit
%%input
a = param.a;
e = param.e;
b = param.b;
mu = G*M;
%%Calculation
%position
r  = R_m+h_atm;
%theta = acos((a*(1-e^2)-r)/(r*e))
%R = r*[cos(theta),sin(theta),0];
x = ((r+a)/-e);
y = sqrt(r^2-x^2);
R = [x,y,0];
%velocity
v = sqrt(mu*(2/r-1/a))
x_rc = (x+a*e)
rc = b^2*x_rc/(a^2*sqrt(b^2*x_rc^2/a^2-b^2))
%rc = b*x_rc/a*sqrt(x_rc^2-a^2)
V_unit = [1,rc,0]/norm([1,rc,0]);
V = v*V_unit;
%acceleration
A = -mu/r^3*R;
%%Output
out.R = R;
out.V = V;
out.A = A;
out
end