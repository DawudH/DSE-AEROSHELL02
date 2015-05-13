function [TriGeom,xvector,yvector,zvector] = Apollo(q)

RShell = 4.694;
Redge = 0.196;

theta = deg2rad(22.8);
thetaedge = pi/2;
thetaedge = linspace(theta,thetaedge+theta,q);
theta = linspace(0,theta,q);
theta = meshgrid(theta);

x = RShell.*cos(theta);
y = RShell.*sin(theta);
x = x(1,:)';
y = y(1,:)';
x2 = 4.138+Redge.*cos(thetaedge);
y2 = 1.956-Redge+Redge.*sin(thetaedge);

x = [x; x2'];
y = [y; y2'];
x = x-RShell;
x = meshgrid(x);
y = meshgrid(y);

phi = linspace(0,2*pi,2*q);
phi = meshgrid(phi);

yr = y.*sin(phi');
zr = y.*cos(phi');

%% Vectorization of grid
xvector = [];
yvector = [];
zvector = [];

for n = 1:length(x(1,:))
    for m = 1:length(x(:,1))
        xvector = [xvector,x(m,n)];
        yvector = [yvector,yr(m,n)];
        zvector = [zvector,zr(m,n)];
    end
end

%% Triangulation matrix 
p = length(xvector);
Tri = [0 0 0];
for i = 0:2*q-2
    for j = 1:2*q-1
%         if i == 0      
%             Tri(2*(i*q+j)-1,:) = [0 0 0];
%             Tri(2*(i*q+j),:) = [i*2*q+j (1+i)*2*q+j+1 (1+i)*2*q+j];
%         else
            Tri(2*(i*(2*q-1)+j)-1,:) = [i*2*q+j (i)*2*q+j+1 (i+1)*2*q+j+1];
            Tri(2*(i*(2*q-1)+j),:) = [i*2*q+j (1+i)*2*q+j+1 (1+i)*2*q+j];
%         end
    end
end
 Tri0 = Tri(:,1) == 0;
 Tri(Tri0,:) = [];
TriGeom = triangulation(Tri, xvector', yvector', zvector');
% FN = faceNormal(TriGeom);
% IC = incenter(TriGeom);
% 
% %trisurf(TriGeom.ConnectivityList,xvector',yvector',zvector');
% hold on
% scatter3(xvector',yvector',zvector')
% axis equal
% quiver3(IC(:,1),IC(:,2),IC(:,3),FN(:,1),FN(:,2),FN(:,3))
