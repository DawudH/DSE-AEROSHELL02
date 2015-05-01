% load constants
constants

ry = 10* R_m; %[m]
v = 7000; %[m/s]
dt = 1;
CD = 1.2;
rx = -4128934.10000000;
CL = 0.25*CD;
tend = 3600 * 24 * 1.5;

[out] = orbit_full(rx,ry,CD,v,dt,CL,tend);