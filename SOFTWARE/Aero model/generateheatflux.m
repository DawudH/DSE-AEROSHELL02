clear;
close all;
clc;
out = open('out.mat');

range = 1:length(out.out.tp);
% range = 1251:4232;

t = out.out.tp;
rho = out.out.rho;
T = out.out.T;
speed_sound = out.out.speed_sound;
M = out.out.M;
V = M.*speed_sound;

gamma = 1.29;
alpha = 10;
q = 24;
[ coords, tri, A ] = generategeometry( 'pastille12m1.5m', q );

qmax = zeros(size(t));
Tboundary = zeros(size(t));
center = [0 0 0];

for i = range
%     i
%     rho(i)
    if rho(i) > 1e-7
        mod = modnewtonian( coords, tri, gamma, speed_sound(i), center, rho(i), T(i), A);
        mod = mod.calcAeroangle(V(i),deg2rad(alpha),0);
        qmax(i) = mod.qmax_array(end);
        Tboundary(i) = mod.Tmax_array(end);
        disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax(i))));  
    end
end

disp('Finished!');
