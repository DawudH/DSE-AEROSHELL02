function dV_dt = diff_V(m,D,g,gamma,r,omega,phi,chi)
dV_dt = -1/m*D-g*sin(gamma)+omega^2*r*cos(phi)^2*(...
    sin(gamma)-cos(gamma)*tan(phi)*sin(chi));
end