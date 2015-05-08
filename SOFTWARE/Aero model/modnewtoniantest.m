
% xvector = [0,0,0,0];
% yvector = [0,1,1,0];
% zvector = [0,0,1,1];
% coords = [xvector;yvector;zvector;];
% mod = modnewtonian( coords, tri, gamma, a, center);


tri = [[1 2 4];[2 3 4]];
a = 300;
gamma = 1.4;
center = zeros(3,1);


q = 5;
r = 3;
R = 10;
t = 3;
type = 't'; %sphere

[TriGeom, xvector, yvector, zvector] = TriMeshGen(q, R, r, t, type);
tri = TriGeom.ConnectivityList;

mod = modnewtonian( TriGeom.Points', tri, gamma, a, center);
mod = mod.calcAeroangle(7e3,deg2rad(10),0);
mod.CRA_array

s = trisurf(tri,xvector,yvector,zvector, mod.Cpdist_array);
axis equal;
colorbar;
xlabel('x')
ylabel('y')
zlabel('z')