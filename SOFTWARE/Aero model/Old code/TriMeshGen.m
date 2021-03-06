function [TriGeom,xvector,yvector,zvector] = TriMeshGen(q,R,r,t,type)

%% Sphere
if type == 's'
%% Define Polar Coordinates
theta = linspace(0,pi,q);              
phi  = linspace(pi,3*pi,q);
[theta,phi]=meshgrid(theta,phi);
%% Shape definition Sphere
x=(r*sin(theta)).*cos(phi)-r;           %x,y,z definitions of a donut in polar. Should be replaced by function??
y=(t*sin(theta)).*sin(phi);
z=R.*cos(theta);
%% Vector setup
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
            Tri(2*(i*q+j)-1,:) = [i*q+j (i+1)*q+j+1 (i)*q+j+1];
            Tri(2*(i*q+j),:) = [0 0 0];
        else
            Tri(2*(i*q+j)-1,:) = [i*q+j (i+1)*q+j+1 (i)*q+j+1];
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j (1+i)*q+j+1];
        end
    end
end
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
id = 'MATLAB:triangulation:PtsNotInTriWarnId';
warning('off',id)
%% Periodic domain; set correct points in triangulaition matrix
for i = 1:q
    Tri(Tri==i*q)=(i-1)*q+1;
end

TriGeom = triangulation(Tri, xvector', yvector', zvector');


%% Torus
elseif type == 't' 
%% Define Polar Coordinates
theta = linspace(-pi,pi,q);              
phi  = linspace(0,2*pi,q);
[theta,phi]=meshgrid(theta,phi);
%% Shape definition Torus
x=r.*cos(theta) -r;
y=(R+t*sin(theta)).*sin(phi);
z=(R+r*sin(theta)).*cos(phi);
%% Vector Setup
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



%% Create triangulation matrix
p = length(xvector);
Tri = [0 0 0];
for i = 0:q-2
    for j = 1:q-1
        Tri(2*(i*q+j)-1,:) = [i*q+j i*q+j+1 (1+i)*q+j+1];
        Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j+1 (1+i)*q+j];
    end
end
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
id = 'MATLAB:triangulation:PtsNotInTriWarnId';
warning('off',id)
TriGeom = triangulation(Tri, xvector', yvector', zvector');    


%% Blunt Half Cone. t = angle of cone, r is half dome radius, R is max radius
elseif type == 'c'
    if r>R 
        print('R must be bigger then r')
    else
X = 1;  
X = linspace(0,X,q);
X = meshgrid(X);

t = 90 -t;
t = tand(t)*r;

theta = linspace(pi,0,q);              
phi  = linspace(pi,3*pi,q);
[theta,phi]=meshgrid(theta,phi);
%% Shape definition Sphere
%x=(r*sin(theta)).*cos(phi);           %x,y,z definitions of a donut in polar. Should be replaced by function??
%y=(r*sin(theta)).*sin(phi);
%z=r.*cos(theta);

%R = 12;
 r = linspace(0,r,q); 
 r = meshgrid(r);
rgrad = max(max(r)); %Gradient of r at the outer edge of the half dome
theta = linspace(0,2*pi,q);
theta = meshgrid(theta);
theta = theta';
 
 %% Cylindrical transformation
 X = (t)*((X.^2)/2);
 z=r.*cos(theta);           
 y=r.*sin(theta);
 x=-X;

%% Coordinate calculation of upper ring(after linear section)
zmax = R.*cos(theta);
ymax = R.*sin(theta);
xmax = -max((t/(R*rgrad))+(R-max(r))*(t/rgrad));
% x5 = -max((t/(2.5*rgrad))+(2.5-max(r))*(t/rgrad))
%xmax = -t/max(max(r))*(R-max(max(r)));
xmax = ones(q)*xmax;
x = [x,xmax(:,1)];
y = [y,ymax(:,1)];
z = [z,zmax(:,1)];

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

xvector(1:q-1)=[];
yvector(1:q-1)=[];
zvector(1:q-1)=[];
p = length(xvector)/q;
for i = 1:p
    xvector(i*q+1) = xvector((i-1)*q+2);
    yvector(i*q+1) = yvector((i-1)*q+2);
    zvector(i*q+1) = zvector((i-1)*q+2);
end


%% Triangulation matrix 
p = length(xvector)/q;
Tri = [0 0 0];
for i = 1:q-1
    for j = 1:q-1
%         if i == 0      
%             Tri(2*(i*q+j)-1,:) = [0 0 0];
%             Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j+1 (1+i)*q+j];
%         else
            Tri(2*(i*q+j)-1,:) = [i*q+j (i)*q+j+1 (i+1)*q+j+1];
            Tri(2*(i*q+j),:) = [i*q+j (1+i)*q+j+1 (1+i)*q+j];
%         end
    end
end
%% Triangulation of base 
Tri1 = [0 0 0];
for I = 2:q
    if I == q
        Tri1(I-1,:) = [1 2 I];
    else
       Tri1(I-1,:) = [1 I+1 I];
    end
end

Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];
Tri = Tri-(q-1);
Tri = [Tri1;Tri];

id = 'MATLAB:triangulation:PtsNotInTriWarnId';
warning('off',id)

for i = 1:q
    Tri(Tri==i*q+1)=(i-1)*q+2;
end

TriGeom = triangulation(Tri, xvector', yvector', zvector');

    end
end





% %% Define geometry in Polar Coordinates
% q = 20;                                 %Number of discrete points used along phi and theta
% R    = 12;                               
% r    = 3;
% theta = linspace(0,2*pi,q);              
% phi  = linspace(0,2*pi,q);
% [theta,phi]=meshgrid(theta,phi);
% 
% % TORUS
% % x=(R+r*cos(theta)).*cos(phi);
% % y=(R+r*cos(theta)).*sin(phi);
% % z=r.*sin(theta);
% 
% % SPHERE
% x=(r*sin(theta)).*cos(phi);           %x,y,z definitions of a donut in polar. Should be replaced by function??
% y=(R*sin(theta)).*sin(phi);
% z=r.*cos(theta);
% 
% 
% %% Create X, Y, Z vectors for grid generations
% xvector = [];
% yvector = [];
% zvector = [];
% for n = 1:length(x(1,:))
%     for m = 1:length(x(:,1))
%         xvector = [xvector,x(m,n)];
%         yvector = [yvector,y(m,n)];
%         zvector = [zvector,z(m,n)];
%     end
% end
% 
% %% For loop to create the triangulation matrix. Tri stores Triangulation matrix entries
% p = length(xvector);
% Tri = [0 0 0];
% for i = 0:q-2
%     for j = 1:q-1
%         Tri(2*(i*q+j)-1,:) = [i*q+j i*q+j+1 q+i*q+j+1];
%         Tri(2*(i*q+j),:) = [i*q+j q+i*q+j+1 q+i*q+j];
%     end
% end
% Tri0 = Tri(:,1) == 0;
% Tri(Tri0,:) = [];
% 
% %% FaceNormal Calculation
% TriGeom = triangulation(Tri, xvector', yvector', zvector');









