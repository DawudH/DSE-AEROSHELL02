clear all
close all

%% Define geometry in Polar Coordinates
q = 100;                                 %Number of discrete points used along phi and theta
R    = 3;                               
r    = 3;
theta= linspace(0,pi,q);              
phi  = linspace(0,2*pi,q);
[theta,phi]=meshgrid(theta,phi);

x=(r*sin(theta)).*cos(phi);           %x,y,z definitions of a donut in polar. Should be replaced by function??
y=(R*sin(theta)).*sin(phi);
z=r.*cos(theta);


%% Create X, Y, Z vectors for grid generations
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

%% For loop to create the triangulation matrix. Tri stores Triangulation matrix entries
p = length(xvector);
Tri = [0 0 0];
for j = 1:(p-(q+1))
    Tri(2*j-1,:) = [j j+1 q+j+1];
    Tri(2*j,:) = [j q+j+1 q+j];
end

%% FaceNormal Calculation
TriGeom = triangulation(Tri, xvector', yvector', zvector');









