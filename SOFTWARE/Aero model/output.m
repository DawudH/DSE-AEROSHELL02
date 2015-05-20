% close all;
clear

shapetexts.horizontalplate = 'horizontalplate';
shapetexts.verticalplate = 'verticalplate';
shapetexts.sphere12m = 'sphere12m';
shapetexts.pastille12m15m = 'pastille12m1.5m';
shapetexts.deg60cone = 'deg60cone';
shapetexts.deg30cone = 'deg30cone';
shapetexts.irvevalidation = 'irvevalidation';
shapetexts.apollovalidation = 'apollovalidation';
shapetexts.torus = 'torus';
shapetexts.concept_irve = 'concept_irve';
shapetexts.concept_apollo = 'concept_apollo';
shapetexts.concept_isotensoid = 'concept_isotensoid';
shapetexts.ballute = 'ballute';

a = 180;
gamma = 1.29;
rho = 9e-5;
T = 130;
q = 37;

alpha0 = -60; %degrees
dalpha = 1; %degrees
alphaend = 60; %degrees
shape = shapetexts.concept_isotensoid;
disp(shape);
disp('Generating geometry...');
[ coords, tri, A, center ] = generategeometry( shape, q );

disp('Initialising...');
mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
disp('Calculating...');
mod = mod.alphasweep(a*20, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

disp('Writing to file...');
dlmwrite('outputfiles/isotensoid.txt', [mod.alpha_array', mod.CRA_aero_array', mod.CMA_aero_array', mod.qmax_array']);
disp('Done!');

