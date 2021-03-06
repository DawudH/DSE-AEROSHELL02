function [TriGeom,xvector,yvector,zvector] = Sharpcone(q,r,hca)
%% Deze code is beun, indien nodig voor later project werk opschonen voor gebruik
%% Define Spherical nose half cone in cylindrical coordinates
% q = 40;
% r = 6;
% rmax = tand(60);

%t=tand(30);  %Gradient of half cone
hca = 90-hca;
R = r*tand(hca);  %Coordinate for lower dome, grid creation
R = linspace(0,R,q);
R = meshgrid(R);

r = linspace(0,r,q); 
r = meshgrid(r);
rgrad = max(max(r)); %Gradient of r at the outer edge of the half dome
theta = linspace(0,2*pi,q);
theta = meshgrid(theta);
theta = theta';

%% Cylindrical transformation
% R = t*(R);
z=r.*cos(theta);           
y=r.*sin(theta);
x=-R;
x = x;
y = y;
z = z;

%% Coordinate calculation of upper ring(after linear section)
% zmax = rmax.*cos(theta);
% ymax = rmax.*sin(theta);
% xmax = -max((t/rgrad)+(rmax-max(r))*(t/rgrad));
% xmax = ones(q)*xmax;
% x = [x,xmax(:,1)];
% y = [y,ymax(:,1)];
% z = [z,zmax(:,1)];

%% Vectorization of grid
xvector = [];
yvector = [];
zvector = [];

for n = 1:length(x(1,:))
    for m = 1:length(x(:,1))
        xvector = [xvector,x(m,n)];
        yvector = [yvector,y(m,n)];
        zvector = [zvector,z(m,n)];
    end
end


%% Triangulation matrix 
p = length(xvector);
Tri = [0 0 0];
for i = 0:q-2
    for j = 1:q-1
        if i == 0      
            Tri(2*(i*q+j)-1,:) = [0 0 0];
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j+1 (1+i)*q+j];
        else
            Tri(2*(i*q+j)-1,:) = [i*q+j (i)*q+j+1 (i+1)*q+j+1];
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j+1 (1+i)*q+j];
        end
    end
end
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
TriGeom = triangulation(Tri, xvector', yvector', zvector');

% trisurf(TriGeom.ConnectivityList,xvector,yvector,zvector)
% axis equal
% xlabel('x')
% ylabel('y')
% zlabel('z')
end







