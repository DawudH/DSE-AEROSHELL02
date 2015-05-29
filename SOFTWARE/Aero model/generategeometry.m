function [ TriGeom, A, center ] = generategeometry(q )
%[ coords, tri, A ] = generategeometry( type, q )
%   Generate used geometries

%     t = 60;
%     R = 6;
%     r = 2;
%     TriGeom = TriMeshGen(q, R, 2, t, 'c');
%     coords = TriGeom.Points';
%     tri = TriGeom.ConnectivityList;
%     A = pi*6^2;
%     t = 90-t;
%     t = tand(t)*r;            
%     center = [-((5/3+t/(2.5*r)+(2.5-r)*t/r)*9000+3/5*(t/R/r+(R-r)*t/r)*1000)/10000, 0, 0];
    %t gradient, r half dome radius, R is max radius            
    
    Nr = q;
    a = 0;
    r = 6;
    h = 3;
    poly = [1 0 0];
    [TriGeom, A] = ParaGeom(Nr, a, r, h, poly);
    center = [0 0 0];
end