% close all;

shapetexts.horizontalplate = 'horizontalplate';
shapetexts.verticalplate = 'verticalplate';
shapetexts.sphere12m = 'sphere12m';
shapetexts.pastille12m15m = 'pastille12m1.5m';
shapetexts.deg60cone = 'deg60cone';
shapetexts.deg30cone = 'deg30cone';
shapetexts.irvevalidation = 'irvevalidation';
shapetexts.apollovalidation = 'apollovalidation';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
q = 75;

alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 10; %degrees

[ coords, tri, A ] = generategeometry( shapetexts.sphere12m, q );

mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
% mod = mod.alphasweep(a*20, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

mod = mod.calcAeroangle(7e3,deg2rad(0),0);

mod.plotCp(true, false);
mod.CRA_aero_array
    
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});

