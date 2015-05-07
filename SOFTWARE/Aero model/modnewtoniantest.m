
x = [0,1,0,0,0,-1];
y = [0,0,1,0,-1,-1];
z = [0,0,0,1,0,1];
coords = [x;y;z;];
tri = [[1,2,3];[1,3,4];[1,4,5];[4,5,6]];
a = 300;
gamma = 1.4;
center = zeros(3,1);

% mod = modnewtonian( coords, tri, gamma, a, center);
mod = modnewtonian( TriGeom.Points', TriGeom.ConnectivityList, gamma, a, center);
mod = mod.calcAeroangle(70000000,deg2rad(45),0);
mod.Cpdist_array;

trisurf(Tri,xvector,yvector, zvector, mod.Cpdist_array);
axis equal;
colorbar;
xlabel('x')
ylabel('y')
zlabel('z')