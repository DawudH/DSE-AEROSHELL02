%no clue why it doesn't work
clc;
clear all;
%number of segments
n = 1000.;
%Mars Radius [km]
R_m = 3390.;
%Mars outline
x_m = linspace(-R_m,R_m,n);
y_m = sqrt(R_m^2 - x_m.^2);
%apogee and perigee height
h_a = 200.;
h_p = 150.;
%semi-major axis, eccentricity and semi-minor axis
x_min = -(R_m+h_p);
x_max = (R_m+h_a);
a = R_m+(h_a+h_p)/2;
e = x_max/a-1;
b = sqrt((1-e^2)*a^2);
x = linspace(x_min,x_max,n);
y = b*sqrt(1-(x/a).^2);
figure
axis normal
plot(x,y,x,-y,x_m,y_m,x_m,-y_m);