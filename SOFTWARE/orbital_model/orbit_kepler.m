function [orbit,t] = orbit_kepler(kepler_param,orbit_init)

a = kepler_param.a;
e = kepler_param.e;
b = kepler_param.b;
theta = kepler_param.theta;
thetap = kepler_param.thetap;
T =  kepler_param.T;

r = norm(orbit_init.R);

% determine time at reentry
    % area of ellipse between exit atmos and going back:
    c = a*e + r * cos(theta)
    Atot = b*a*pi;
    A = 1/2*b*a*pi + b*c*sqrt(1-c^2/a^2) - b*a*asin(-c/a) - 1/2*r^2*cos(theta)*sin(theta);
    
    t = A/Atot * T;

% Determine the R, V, and accelerations at the end point:
    angle = 2*theta;
    Tz = [      cos(angle),     sin(angle), 0;...
                -sin(angle),    cos(angle), 0;...
                0,              0,          1];
    
    % rotation and translation matrix for V
    
        % convert to new point..
        Ve = Tz * orbit_init.V';
        % rotate to another axis system:
        angle = thetap;
        Tv = [  cos(angle),     sin(angle), 0;...
                -sin(angle),    cos(angle), 0;...
                0,              0,          1];
        Ve = Tv * Ve;     
        %reflect around ye axis
        Refy = [1, 0, 0;...
                 0, -1, 0;...
                 0, 0, 1];
        Ve = Refy * Ve;
        % rotate back to express in 0 frame
        Tv = [  cos(angle),     -sin(angle), 0;...
                sin(angle),    cos(angle), 0;...
                0,              0,          1];
        Ve = Tv * Ve;
    R = (Tz * orbit_init.R')'
    V = Ve'
    acc = (Tz * orbit_init.a')'
    ag = (Tz * orbit_init.ag')'

    

% output        
orbit.R = R;
orbit.V = V;
orbit.a = acc;
orbit.ad = [0,0,0];
orbit.al = [0,0,0];
orbit.ag = ag;

end