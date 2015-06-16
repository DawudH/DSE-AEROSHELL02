clear;
close all;
clc;
casename = 'iteration1_1_orbit';
load(strcat('aeroshapes\',casename,'.mat'));
skewnessoriginal = x(1);
heightfactor = x(2);

poly = [x(3:end),0,0];   

for diameter = 18
    diameter
    load(strcat('orbits\out_d', num2str(diameter), '_one_time.mat'));    
    clear mod;
    radius = diameter/2;
    height = radius*heightfactor;


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
    qmax_array = zeros(size(t));
    Tboundary = zeros(size(t));

    [ mod, center ] = generateGeometry(poly, q, skewness, radius, height, speed_sound(1), gamma, rho(1), T(1));
    geom = mod.geom;
    

    parfor i = range
        if rho(i) > 1e-14
            mod = modnewtonian(geom, gamma, speed_sound(i), center, rho(i), T(i));
            mod = mod.calcAeroangle(V(i),deg2rad(alpha(i)), beta, phi);
            [Tmax, qmax] = mod.calcStagnationHeatFlux();
            qmax_array(i) = qmax;
            Tboundary(i) = Tmax;
%             disp(strcat('current number: (', num2str(i), '/', num2str(max(range)), '), qmax: ', num2str(qmax_array(i))));  
        end
    end
    A_whetted = sum(geom.areas);
    disp('Finished!');
    save(strcat('heatflux\heatflux_out_d', num2str(diameter), '_one_time.mat'));
    plot(t, qmax_array);
end