clear;
close all;
clc;
casename = 'orbit_iteration_1_1_entry';
load(strcat('orbits\',casename,'.mat'));
load('aeroshapes\iteration1_1_heatflux.mat');

range = 1:length(out.tp);

t = out.tp;
rho = out.rho;
T = out.T;
speed_sound = out.speed_sound;
M = out.M;
V = M.*speed_sound;

gamma = 1.29;
alpha = out.alpha;
beta = 0;
phi = 0;
q = 31;
radius = 15;
qmax_array = zeros(size(t));
Tboundary = zeros(size(t));

Nr = q;
a = 0;
h = 0.5*radius;
poly = [1 0 0];
center = [0 0 0];

% [ TriGeom, A, center ] = generategeometry( q, radius );
% geom = aeroGeometry(TriGeom, A, poly);
% mod = modnewtonian(geom, gamma, a, center, rho, T);
geom = mod.geom;

parfor i = range
    if rho(i) > 1e-14
        mod = modnewtonian(geom, gamma, speed_sound(i), center, rho(i), T(i));
        mod = mod.calcAeroangle(V(i),deg2rad(alpha(i)), beta, phi);
        [Tmax, qmax] = mod.calcStagnationHeatFlux();
        Tmax = Tmax(end);
        qmax = qmax(end);
        qmax_array(i) = qmax;
        Tboundary(i) = Tmax;
        disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax_array(i))));  
    end
end
A_whetted = sum(geom.areas);
disp('Finished!');
save(strcat('heatflux/heatflux_',casename,'.mat'));
plot(t, qmax_array);