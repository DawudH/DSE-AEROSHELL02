
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

a = 300;
gamma = 1.4;
rho = 1e-5;
T = 150;
q = 31;

alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 40; %degrees

[ coords, tri, A, center ] = generategeometry( shapetexts.verticalplate, q );

mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
% mod = mod.alphasweep(a*20, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

mod = mod.calcAeroangle(7e3,deg2rad(30),deg2rad(0));

mod.plotCp(true, false);
% mod.CR_aero_array
    
  
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});

