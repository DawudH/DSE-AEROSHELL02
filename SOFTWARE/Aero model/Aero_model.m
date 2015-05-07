V_inf=[0,0,3000];
R    = 5;
r    = 3;

theta= linspace(0,2*pi,10);
phi  = linspace(0,2*pi,10);


[theta,phi]=meshgrid(theta,phi);
x=(R+r*cos(theta)).*cos(phi);
y=(R+r*cos(theta)).*sin(phi);
z=r.*sin(theta);
kleurtjes = z;

surface=surf(x,y,z, kleurtjes);

axis equal
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



% dt = DelaunayTri(xvector', yvector', zvector');
% [tri,Xb] = freeBoundary(dt);
% figure;
% trisurf(tri,Xb(:,1),Xb(:,2),Xb(:,3), 'FaceColor', 'cyan', 'faceAlpha', 0.8);

