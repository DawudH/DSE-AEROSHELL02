% load('aeroshapes\iteration1_0_simplified.mat');
clearvars -except x


q = 61;


params = globalParams();

skewness = x(1);

heightfactor = x(2);
height = params.radius * heightfactor;


poly = [x(3:end),0,0];

alpha0 = -40;
dalpha = 1;
alphaend = 40;
V = 7000;
beta = 0;
phi = 0;

[mod, center] = generateGeometry(poly, q, skewness, params.radius, height);
mod = mod.alphasweep(V, beta, phi, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));
save('aeroshapes\iteration1_0_heatflux.mat');
