close all;
% xvector = [0,0,0,0,-1,-1];
% yvector = [0,1,1,0,1,0];
% zvector = [0,0,1,1,0,0];
% coords = [xvector;yvector;zvector;];


tri = [1 2 5; 1 5 6];%;1 2 4;2 3 4;1 5 2;1 6 5];
a = 200;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;

% mod = modnewtonian( coords, tri, gamma, a, center);


q = 100;
r = 3;
R = 12;
t = 12;
type = 's'; %sphere

[TriGeom, xvector, yvector, zvector] = TriMeshGen(q, R, r, t, type);
tri = TriGeom.ConnectivityList;

mod = modnewtonian( TriGeom.Points', tri, gamma, a, center, rho, T);
mod = mod.alphasweep(7e3, 0, deg2rad(0), deg2rad(20), deg2rad(2));
% mod = mod.calcAeroangle(7e3,deg2rad(10),0);
mod.plotCp(true, true);
mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});