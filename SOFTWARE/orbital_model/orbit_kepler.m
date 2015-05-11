function [orbit,t] = orbit_kepler(kepler_param,orbit_init)

constants

a = kepler_param.a;
e = kepler_param.e;
b = kepler_param.b;
theta = kepler_param.theta;
thetap = kepler_param.thetap;

T =  kepler_param.T;

r = norm(orbit_init.R);
if a>0
    % determine time at reentry
        % area of ellipse between exit atmos and going back:
        c = a*e + r * cos(theta);
        Atot = b*a*pi;
        A = 1/2*b*a*pi + b*c*sqrt(1-c^2/a^2) - b*a*asin(-c/a) - 1/2*r^2*cos(theta)*sin(theta);

        t = A/Atot * T;

        % convert V to axis system with R as X axis..
            angle = thetap;
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
            acc = Te' * ae;
else
    r1 = R_m+h_atm-1
    theta1 = acos((a*(1-e^2)-r1)/(r1*e))
    theta2 = theta1+2*pi/3600
    r2 = a*(e^2-1)/(1+e*cos(theta))
    R = r1*[cos(theta1),sin(theta1),0]
    R2 = r2*[cos(theta2),sin(theta2),0]
    v = sqrt(G*M_mars*(2/r1-1/a))
    V = (R2-R)/norm(R2-R)*v
    acc = -G*M_mars/r^3*R
    t=0;
end

% output        
orbit.R = R;
orbit.V = V;
orbit.a = acc;
orbit.ad = [0,0,0];
orbit.al = [0,0,0];
orbit.q = 0;
orbit.speed_sound = 0;
orbit.M = 0;
orbit.ag = acc;

end