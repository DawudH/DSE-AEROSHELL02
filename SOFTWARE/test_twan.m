clear all
clc

M_mars = 6.419*10^23; %[kg]
G = 6.673*10^-11; %[N*(m/kg)^2]
R_m = 6794000/2; %[m]
x_m_h = linspace(-R_m,R_m,1000);
x_m_r = flip(x_m_h);
x_m = [x_m_h,x_m_r];
y_m_p = sqrt(R_m^2-x_m_h.^2);
y_m = [y_m_p,-y_m_p];
h = 150000; %[m]
r = R_m+h; %[m]
x_a_h = linspace(-r,r,1000);
x_a_r = flip(x_a_h);
x_a = [x_a_h,x_a_r];
y_a_p = sqrt(r^2-x_a_h.^2);
y_a = [y_a_p,-y_a_p];
v = 6030; %[m/s]
%v = sqrt(G*M_mars*2/r);
a = (2/r-v^2/(G*M_mars))^-1;
b = 3000000;
if a<0
    x_h = linspace(-a,a,1000);
    x_h_r = flip(x_h);
    x = [x_h,x_h_r];
    y_p = b*sqrt(((x_h-2*a)/a).^2-1);
    y_n = -b*sqrt(((x_h_r-2*a)/a).^2-1);
    y = [y_p,y_n];
elseif a>0
    x_h = linspace(-a,a,1000);
    x_h_r = flip(x_h);
    x = [x_h,x_h_r];
    y_p = b*sqrt(1-(x_h/a).^2);
    y_n = -b*sqrt(1-(x_h_r/a).^2);
    y = [y_p,y_n];
else
    fprintf('parabola');
end
plot(x,y,x_m,y_m,x_a,y_a);
