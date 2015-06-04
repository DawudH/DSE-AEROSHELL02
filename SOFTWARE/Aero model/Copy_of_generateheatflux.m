clear;
close all;
clc;
casename = 'out_d6_one_time';
load(strcat('orbits\',casename,'.mat'));

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
radius = 3;
qmax_array = zeros(size(t));
Tboundary = zeros(size(t));
a = 300;

[ TriGeom, A, center ] = generategeometry( q, radius );
geom = aeroGeometry(TriGeom, A);
mod = modnewtonian(geom, gamma, a, center, rho, T);

parfor i = range
    if rho(i) > 1e-14
        mod = modnewtonian(geom, gamma, speed_sound(i), center, rho(i), T(i));
        mod = mod.calcAeroangle(V(i),deg2rad(alpha(i)), beta, phi);
        [Tmax, qmax] = mod.calcStagnationHeatFlux();
        qmax_array(i) = qmax;
        Tboundary(i) = Tmax;
        disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax_array(i))));  
    end
end

disp('Finished!');
save(strcat('heatflux/heatflux_',casename,'.mat'));