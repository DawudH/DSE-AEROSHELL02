clear;
close all;
clc;
out = open('orbits/orbit.mat');

range = 1:length(out.out.tp);

t = out.out.tp;
rho = out.out.rho;
T = out.out.T;
speed_sound = out.out.speed_sound;
M = out.out.M;
V = M.*speed_sound;

gamma = 1.29;
alpha = 30;
q = 22;
% [ coords, tri, A ] = generategeometry( shapetexts.concept_isotensoid, q );

qmax_array = zeros(size(t));
Tboundary = zeros(size(t));
center = [0 0 0];
a = 300;

[ TriGeom, A, center ] = generategeometry( q );
geom = aeroGeometry(TriGeom, A);
mod = modnewtonian(geom, gamma, a, center, rho, T);
Tw = zeros(size(mod.geom.areas));

parfor i = range

    if rho(i) > 1e-10
%         [ TriGeom, A, center ] = generategeometry( q );
%         geom = aeroGeometry(TriGeom, A);
        mod = modnewtonian(geom, gamma, speed_sound(i), center, rho(i), T(i));
%         mod.a_inf = speed_sound(i);
%         mod.rho_inf = rho(i);
%         mod.T_inf = T(i);
        mod = mod.calcAeroangle(V(i),deg2rad(alpha),0);
        [Tmax, qmax, qw] = mod.calcStagnationHeatFlux(Tw);
        qmax_array(i) = qmax;
        Tboundary(i) = Tmax;
        disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax_array(i))));  
    end
end

disp('Finished!');
