clear all
clc
comment = false;

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
h_a = 150000;
r_a = R_m+h_a;
x_a_h = linspace(-r_a,r_a,1000);
x_a_r = flip(x_a_h);
x_a = [x_a_h,x_a_r];
y_a_p = sqrt(r_a^2-x_a_h.^2);
y_a = [y_a_p,-y_a_p];
v = 7000; %[m/s] CHANGE
%v = sqrt(G*M_mars*2/r);
a = (-v^2/(G*M_mars))^-1;
e = 1.35; %CHANGE
b = sqrt(a^2*(e^2-1));
if comment
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
    figure;
    plot(x,y,x_m,y_m,x_a,y_a);
end
h = 2000000;
r = R_m + h;
theta = 40/180*pi;
dia = 12;
S = pi/4*dia^2;
m = 10000; %[kg]
CD = 1; %[-]
r_v = [-r*sin(theta) r*cos(theta) 0];
r_v_n = r_v;
v_v = [0 -v 0];
dt = 1/10;
orbit = true;
t=0;
while orbit == true
    if (abs(r)-R_m)<150000
        rho = 1;
    else
        rho = 0;
    end
    a_v = -(G*M_mars/r^3)*r_v_n - CD*0.5*rho*v_v.^2*S/m;
    r_v_n = r_v_n + v_v * dt + a_v *dt^2;
    r_v = [r_v r_v_n];
    v_v = v_v + a_v * dt;
    t = t+dt;
    if (abs(r_v)<R_m) | (abs(r_v) > 2*R_m)
        orbit = false;
    end
end
x = r_v(1 : 3 : end);
y = r_v(2 : 3 : end);
z = r_v(3 : 3 : end);
[x_m, y_m, z_m] = sphere;
x_m = R_m * x_m;
y_m = R_m * y_m;
z_m = R_m * z_m;
figure;
hold on
surf(x_m,y_m,z_m);
scatter3(x,y,z);