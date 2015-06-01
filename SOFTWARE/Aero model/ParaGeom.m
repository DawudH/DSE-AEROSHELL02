function [TriGeom, AFrontal] = ParaGeom(Nr,a,r,h,poly)

AFrontal = pi*r^2;
% Nr = 40;
%% Polar coordinate definition
theta = linspace (0,2*pi,Nr);
r = linspace(0,r,Nr);
c = linspace(0,1,Nr);
% a = 0.5;
x = zeros(Nr);
y = zeros(Nr);
z = zeros(Nr);

%% Geometry creation. 
% Generates a set of circles. If a == 0, these will be concentric around the
% origin. If not, the circle centers will shift away from the origin
% following aa distribution defined by c. This will introduce a skew in the
% shape.
% The circle height above the xy-plane is then determined by z. Adjusting
% this distribution will alter the unskewed shape of the aeroshell. 

for i = 1:Nr
z(i,:) = r(i).*cos(theta)-c(i)*a;
y(i,:) = r(i).*sin(theta);
zx = linspace((i-1)/Nr,(i-1)/Nr,Nr);
x(i,:) = -h*polyval(poly,zx)/sum(poly);
end


%% Vectorization of the coordinates
%this section vectorizes the coordinate matrices. The duplicate points in the origin
%are also removed. 

xvector = [];
yvector = [];
zvector = [];
for n = 1:length(x(1,:))
    for m = 1:length(x(:,1))
        xvector = [xvector,x(n,m)];
        yvector = [yvector,y(n,m)];
        zvector = [zvector,z(n,m)];
    end
end

xvector = xvector';
yvector = yvector';
zvector = zvector';

xvector(1:Nr-1) = [];
yvector(1:Nr-1) = [];
zvector(1:Nr-1) = [];

%% Triangulation Matrix
% Matrix that relates the various points that need to be connected by the
% triangular mesh to each other.

%General triangulation
Tri = [0 0 0];
for j = 0:Nr-3
    for i = 1:Nr-1
        Tri(2*(i+j*Nr)-1,:) = [j*Nr+i+1 j*Nr+i+2 (j+1)*Nr+i+2];
        Tri(2*(i+j*Nr),:) = [j*Nr+i+1 (j+1)*Nr+i+2 (j+1)*Nr+i+1];
    end
end
Tri0 = Tri(:,1) == 0;
Tri(Tri0,:) = [];

% Replace the periodic points by their exact counterparts.
for i = 1:Nr-1
   TriDup2 = Tri(:,2) == i*Nr+1;
   TriDup3 = Tri(:,3) == i*Nr+1;
   Tri(TriDup2,2) = (i-1)*Nr+2;
   Tri(TriDup3,3) = (i-1)*Nr+2;
   xvector(i*Nr+1) = 0 ;
   yvector(i*Nr+1) = 0 ;
   zvector(i*Nr+1) = 0 ;
end
IndDisc = find(~xvector);
IndDisc(IndDisc==1)=[];
xvector(IndDisc) = [];
yvector(IndDisc) = [];
zvector(IndDisc) = [];
Tri1 = Tri(:,1);
Tri2 = Tri(:,2);
Tri3 = Tri(:,3);
for i = 1:length(IndDisc)
    Tri1(Tri1>IndDisc(i)) = Tri1(Tri1>IndDisc(i))-1;
    Tri2(Tri2>IndDisc(i)) = Tri2(Tri2>IndDisc(i))-1;
    Tri3(Tri3>IndDisc(i)) = Tri3(Tri3>IndDisc(i))-1;
    IndDisc = IndDisc -1;
end

Tri = [Tri1 Tri2 Tri3];
% Create triangulation around the origin point. 
TriOrigin = [0 0 0];
for I = 2:Nr
    if I == Nr
        TriOrigin(I-1,:) = [1 2 I];
    else
       TriOrigin(I-1,:) = [1 I+1 I];
    end
end
Tri = [TriOrigin; Tri];
TriGeom = triangulation(Tri, xvector, yvector, zvector);

    









