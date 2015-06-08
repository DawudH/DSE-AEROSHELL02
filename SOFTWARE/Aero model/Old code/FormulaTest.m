clear all 
close all

%% Define Polar Coordinates
q=5; 

t = linspace(0,1,q);
r = linspace(0,1,q);

[r,t] = meshgrid(r,t);

theta = linspace(0,pi,q);              
phi  = linspace(0,pi,q);
[theta,phi]=meshgrid(theta,phi);

R = meshgrid(linspace(0,1,q));
R = (0.2)*R.^2;
%% Shape definition Cylinder
x=r.*cos(phi);           %x,y,z definitions of a donut in polar. Should be replaced by function??
y=r.*sin(phi);
z=R;%.*cos(theta) ;
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
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j (1+i)*q+j+1];
        elseif i == q-2
            Tri(2*(i*q+j)-1,:) = [i*q+j (1+i)*q+j+1 i*q+j+1];
            Tri(2*(i*q+j),:) = [0 0 0];
        else
            Tri(2*(i*q+j)-1,:) = [i*q+j (1+i)*q+j+1 i*q+j+1];
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j (1+i)*q+j+1];
        end
    end
end
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
TriGeom = triangulation(Tri, xvector', yvector', zvector');


trisurf(TriGeom.ConnectivityList,xvector,yvector,zvector)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')