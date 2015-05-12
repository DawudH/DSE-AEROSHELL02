% close all;
horizontalplate = 'horizontalplate';
verticalplate = 'verticalplate';
sphere12m = 'sphere12m';
pastille12m15m = 'pastille12m1.5m';
deg60cone = 'deg60cone';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
q = 70;

alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 15; %degrees

[ coords, tri, A ] = generategeometry( deg60cone, q );

mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
% mod = mod.alphasweep(7e3, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

mod = mod.calcAeroangle(7e3,deg2rad(0),0);

mod.plotCp(true, true);
mod.CR_aero_array
    
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});

