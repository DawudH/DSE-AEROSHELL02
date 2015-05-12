function [ out ] = ek_Parameters(R_init,V_init,A_init,G,M,dt)
%Calculates the kepler orbit parameters needed to calculate the eliptic
%flightpath

    %%Input
    %initial position
    r_init = norm(R_init);
    %position after dt
    R1 = R_init+V_init*dt+A_init*dt^2;
    r1 = norm(R1);
    %initial velocity
    v_init = norm(V_init);
    %angle change between R_init and R1
    theta0 = atan2(R_init(2),R_init(1));
    theta1 = atan2(R1(2),R1(1));
    %gravitational constant of Mars
    mu = G*M;

    %%Calculations
    %determine semi-major axis
    a = (2/r_init-v_init^2/mu)^-1;
    % determine e with areal velocity
    dtheta = theta1 - theta0;
    dA = 1/2 * r_init * r1 * sin(dtheta);
    e = sqrt( 1 - 4 * (dA / dt)^2 / (a * mu) );
    %Determine semi-minor axis
    b = a*sqrt(1-e^2);
    %Perifocal distance
    rp  = a*(1-e);
    %Apofocal distance
    ra = 2*a - rp;
    %Calculate theta (wrt the elipse reference frame [Xe,Ye,Ze])
    theta = acos( (a*(1-e^2)-r_init) / (e*r_init) );
    %Calculate the angle between X0 and Xe
    theta_p = theta0 - theta;
    % calculate the orbital period
    T = 2*pi * sqrt(a^3 / mu);

    %%Output
    out.a = a;
    out.e = e;
    out.b = b;
    out.rp = rp;
    out.ra = ra;
    out.theta = theta;
    out.theta0 = theta0;
    out.theta_p = theta_p;
    out.T = T;
end