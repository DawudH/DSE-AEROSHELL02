function [out] = initial(r,v,theta0,gamma,mu)

out.R = r*[cos(theta0),sin(theta0),0];
out.V = (rotz(gamma)*v*[-sin(theta0),cos(theta0),0]')';
out.A = -mu/r^3*out.R;
out.speed_sound = 0;
out.M = 0;
out.Ag = out.A;
out.Ad = [0,0,0];
out.Al = [0,0,0];
out.J = [0,0,0];
out.q = 0;
out.T = 4; %[K]
out.rho = 0;