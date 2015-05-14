function [ out ] = ek_Endpoint(param,orbit_init)
%Calculates the endpoint of the hyperbolic entry orbit
    
    %%Input
    a = param.a;
    e = param.e;
    theta = param.theta;
    theta_p = param.theta_p;
    r = norm(orbit_init.R);
    b = param.b;
    T = param.T;
    c = a*e + r * cos(theta);
    Atot = b*a*pi;
    A = 1/2*b*a*pi + b*c*sqrt(1-c^2/a^2) - b*a*asin(-c/a) - 1/2*r^2*cos(theta)*sin(theta);

    out.t_kep = A/Atot * T;

    % convert V to axis system with R as X axis..
    angle = theta_p;
    Te = [      cos(angle),     sin(angle), 0;...
                -sin(angle),    cos(angle), 0;...
                0,              0,          1];
    Re = Te * orbit_init.R';
    Ve = Te * orbit_init.V';
    ae = Te * orbit_init.a';
    %reflect around yr axis
    Refy = [ 1,   0, 0;...
             0,   -1,0;...
             0,   0, 1];
    Re = Refy * Re;
    Ve = Refy * Ve;
    ae = Refy * ae;
    % convert back to 0-frame
    R = Te' * Re;
    V = -(Te' * Ve); %change direction of V
    A = Te' * ae;

    %%Output
    out.R = R';
    out.speed_sound = 0;
    out.V = V';
    out.M = 0;
    out.A = A';
    out.Ag = A';
    out.Ad = [0,0,0];
    out.Al = [0,0,0];
    out.J = [0,0,0];
    out.q = 0;
    out.T = 0;
    out.rho = 0;
end
