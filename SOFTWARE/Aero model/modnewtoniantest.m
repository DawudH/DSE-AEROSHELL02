clear;
a = 300;
gamma = 1.4;
rho = 1e-5;
T = 150;
q = 10;

alpha0 = -25; %degrees
dalpha = 5; %degrees
alphaend = 25; %degrees

[ TriGeom, A, center ] = generategeometry( q );
geom = aeroGeometry(TriGeom, A);
% tic
% geom.calcDistances(geom.coords);
% toc
mod = modnewtonian(geom, gamma, a, center, rho, T);
% mod = mod.alphasweep(a*20, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

mod = mod.calcAeroangle(2000,deg2rad(20),deg2rad(0));

Tw = 500*ones(size(mod.Cpdist_array(:,end)));
[Tmax, qmax, qw] = mod.calcStagnationHeatFlux(Tw);

geom.plotValues(qw, 'qw', [0 12], true, false);

% mod.plotCp(true, false);
% mod.CR_aero_array
    
  
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});

