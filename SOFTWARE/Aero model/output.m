% close all;

shapetexts.horizontalplate = 'horizontalplate';
shapetexts.verticalplate = 'verticalplate';
shapetexts.sphere12m = 'sphere12m';
shapetexts.pastille12m15m = 'pastille12m1.5m';
shapetexts.deg60cone = 'deg60cone';
shapetexts.deg30cone = 'deg30cone';
shapetexts.irvevalidation = 'irvevalidation';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
q = 80;

alpha0 = -60; %degrees
dalpha = 1; %degrees
alphaend = 60; %degrees

[ coords, tri, A ] = generategeometry( shapetexts.pastille12m15m, q );

mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
mod = mod.alphasweep(a*20, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

dlmwrite('outputfiles/pastille.txt', [mod.alpha_array', mod.CRA_aero_array', mod.CMA_aero_array']);


