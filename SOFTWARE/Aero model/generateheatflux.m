clear;
close all;
clc;
out = open('out.mat');

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

range = 1:length(out.out.tp);
% range = 1251:4232;

t = out.out.tp;
rho = out.out.rho;
T = out.out.T;
speed_sound = out.out.speed_sound;
M = out.out.M;
V = M.*speed_sound;

gamma = 1.29;
alpha = 30;
q = 7;
[ coords, tri, A ] = generategeometry( shapetexts.concept_irve, q );

qmax = zeros(size(t));
Tboundary = zeros(size(t));
center = [0 0 0];
mod = modnewtonian( coords, tri, gamma, speed_sound(1), center, rho(1), T(1), A);
for i = range
%     i
%     rho(i)
    if rho(i) > 1e-7
%         mod = modnewtonian( coords, tri, gamma, speed_sound(i), center, rho(i), T(i), A);
        mod.a_inf = speed_sound(i);
        mod.rho_inf = rho(i);
        mod.T_inf = T(i);
        mod = mod.calcAeroangle(V(i),deg2rad(alpha),0);
        qmax(i) = mod.qmax_array(end);
        Tboundary(i) = mod.Tmax_array(end);
        disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax(i))));  
    end
end

disp('Finished!');
