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
R = r*[cos(theta),sin(theta),0];
%Express in 0-reference frame
R0 = rotz(rad2deg(theta_p))*R';
x = R0(1);

%Velocity
%magnitude
v = sqrt(mu*(2/r-1/a));
%direction
x_rc = (x+a*e);
rc = a*x_rc/(b^2*sqrt(x_rc^2/b^2+1));
V_unit = [-rc,1,0]/norm([-rc,1,0]);
V = v*V_unit;

%Acceleration
A = -mu/r^3*R;

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

%Plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
figure('name','Orbit')
grid on
axis equal
hold on
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g');
theta_plot = param.theta:0.001:theta;
rk = a * (1- e^2) ./ (1 + e * cos(theta_plot));
polar(theta_plot+param.theta_p,rk,'k');
plot(param.rp*cos(param.theta_p),param.rp*sin(param.theta_p),'*');
plot(-param.ra*cos(param.theta_p),-param.ra*sin(param.theta_p),'d');
end