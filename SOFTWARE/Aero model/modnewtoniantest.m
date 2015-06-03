clear;


a = 150; %speed of sound
gamma = 1.29; %gamma
rho = 1e-5; %Density
T = 150; %Temperatuure
q = 80; % Maat voor aantal elementen
M = 40; %Mach number
V = a*M; %Velocity
beta = 0; %Sideslip
alpha = 0; %Angle of attack

[ TriGeom, A, center ] = generategeometry( q );
geom = aeroGeometry(TriGeom, A);

mod = modnewtonian(geom, gamma, a, center, rho, T);
mod = mod.calcAeroangle(V, alpha, beta);

plot(mod.alpha_array, mod.CRA_aero_array)





% Tw = 500*ones(size(mod.Cpdist_array(:,end)));
% [Tmax, qmax, qw] = mod.calcStagnationHeatFlux(Tw);

% geom.plotValues(qw, 'qw', [0 max(qw)], true, false);

% mod.plotCp(true, false);
% mod.CR_aero_array
    
  
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});
