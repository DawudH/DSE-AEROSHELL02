V_inf=[0,0,3000];
R    = 5;
r    = 3;
theta= linspace(0,2*pi,100);
phi  = linspace(0,2*pi,100);

[theta,phi]=meshgrid(theta,phi);
x=(R+r*cos(theta)).*cos(phi);
y=(R+r*cos(theta)).*sin(phi);
z=r.*sin(theta);
kleurtjes = z;

surface=surf(x,y,z, kleurtjes);

axis equal
