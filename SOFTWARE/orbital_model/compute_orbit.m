% load constants
constants

ry = 10* R_m; %[m]
v = 7000; %[m/s]
dt = 1;
CD = 0.8;
rx = -4133000;
CL = 0.2*CD;
tend = 3600 * 24 * 10;

[out] = orbit_full(rx,ry,CD,v,dt,CL,tend);