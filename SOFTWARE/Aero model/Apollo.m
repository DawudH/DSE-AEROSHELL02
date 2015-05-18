function [TriGeom,xvector,yvector,zvector] = Apollo(q)
% close all
% clear all
% q=7;
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

phi = linspace(0,pi,q);
phi2 = linspace(pi,2*pi,q);
phi = [phi phi2];
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

xvector(1:2*q-1)=[];
yvector(1:2*q-1)=[];
zvector(1:2*q-1)=[];
% xvector(q+1)=[];
% yvector(q+1)=[];
% zvector(q+1)=[];
p = length(xvector)/(q*2);



for i = 1:p
    xvector(i*2*q+1) = xvector((i-1)*2*q+2);
    yvector(i*2*q+1) = yvector((i-1)*2*q+2);
    zvector(i*2*q+1) = zvector((i-1)*2*q+2);
end

Vector = [xvector' yvector' zvector'];

for i = 1:length(Vector)-1
    if Vector(i,:) == Vector (i+1,:)
       Vector(i+1,:) = [0 0 0]; 
    end
end

VectorI = Vector(:,1) == 0 & Vector(:,2) == 0 & Vector(:,2) == 0;
Vector(VectorI,:)=[];
xvector = Vector(:,1);
yvector = Vector(:,2);
zvector = Vector(:,3);

%% Triangulation matrix 

Tri = [0 0 0];
for i = 1:2*q-2
    for j = 1:2*q-1
%         if i == 0      
%             Tri(2*(i*(2*q-1)+j)-1,:) = [0 0 0];
%             %Tri(2*(i*(2*q-1)+j),:) = [i*2*q+j (1+i)*2*q+j+1 (1+i)*2*q+j];
%         else
            Tri(2*(i*(2*q-1)+j)-1,:) = [i*2*q+j (i)*2*q+j+1 (i+1)*2*q+j+1];
            Tri(2*(i*(2*q-1)+j),:) = [i*2*q+j (1+i)*2*q+j+1 (1+i)*2*q+j];
%          end
    end
end
Tri1 = [0 0 0];
for I = 2:2*q
    if I ==2*q
        Tri1(I-1,:) = [1 2 I];
    else
        Tri1(I-1,:) = [1 I+1 I];
    end
end
Tri1(q,:)=[];
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
Tri = Tri-(2*q-1);
Tri = [Tri1;Tri];

 
 
TriGeom = triangulation(Tri, xvector, yvector, zvector);
% FN = faceNormal(TriGeom);
% IC = incenter(TriGeom);
% 
% trisurf(TriGeom.ConnectivityList,xvector',yvector',zvector');
% hold on
% scatter3(xvector',yvector',zvector')
% axis equal
% quiver3(IC(:,1),IC(:,2),IC(:,3),FN(:,1),FN(:,2),FN(:,3))
