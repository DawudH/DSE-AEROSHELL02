function [ coords, tri, A, center ] = generategeometry( type, q )
%[ coords, tri, A ] = generategeometry( type, q )
%   Generate used geometries


    switch type
        case 'torus'
            R = 6;
            r = R/8;
            TriGeom = TriMeshGen(q, R, r, r, 't');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*((R+r)^2-(R-r)^2);
            center = [0 0 0];
        case 'horizontalplate'
            xvector = [0 1 0 1];
            yvector = [0 0 1 1];
            zvector = [0 0 0 0];
            coords = [xvector;yvector;zvector];
            tri = [1 2 3;2 4 3];
            A = 1;
            center = [0 0 0];
        case 'verticalplate'
            xvector = [0 0 0 0];
            yvector = [0 0 1 1];
            zvector = [0 1 1 0];
            coords = [xvector;yvector;zvector];        
            tri = [1 4 2;2 4 3];
            A = 1;
            center = [0 0 0];
        case 'sphere12m'
            TriGeom = TriMeshGen(q, 6, 6, 6, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;
            center = [0 0 0];
        case 'concept_isotensoid'
            r = 5;
            R = 6;
            TriGeom = TriMeshGen(q, R, r, R, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;   
            center = [(r/2*1000+(5/3*+5)*9000)/10000, 0, 0];
        case 'pastille12m1.5m'
            TriGeom = TriMeshGen(q, 6, 1.5, 6, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList; 
            A = pi*6^2;
            center = [0 0 0];
        case 'concept_apollo'
            TriGeom = TriMeshGen(q, 2.5, 0.564960922454593, 2.5, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList; 
            A = pi*2.5^2;  
            center = [-1.67 0 0];
        case 'deg60cone'
            [TriGeom,xvector,yvector,zvector] = Sharpcone(q,6,30);
            coords = [xvector;yvector;zvector];
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;
            center = [0 0 0];
        case 'deg30cone'
            [TriGeom,xvector,yvector,zvector] = Sharpcone(q,6,15);
            coords = [xvector;yvector;zvector];
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;
            center = [0 0 0];
        case 'irvevalidation'
            TriGeom = TriMeshGen(q, 2.93/2, 0.21, 60, 'c');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*1/4*2.93^2;
            center = [0 0 0];
            %t gradient, r half dome radius, R is max radius
        case 'concept_irve'
            TriGeom = TriMeshGen(q, 6, 2, 60, 'c');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;
            center = [0 0 0];
            %t gradient, r half dome radius, R is max radius            
        case 'apollovalidation'
            [TriGeom, xvector, yvector, zvector] = Apollo(q);
            coords = [xvector;yvector;zvector];
            tri = TriGeom.ConnectivityList;
            A = pi*4.694^2;
            center = [0 0 0];
        otherwise
            warning(strcat('The following type is not supported: ', num2str(type)));
            xvector = [0,0,0,0,-1,-1];
            yvector = [0,1,1,0,1,0];
            zvector = [0,0,1,1,0,0];
            coords = [xvector;yvector;zvector;];
            tri = [1 2 4;2 3 4;1 5 2;1 6 5];
            A = 1;
            center = [0 0 0];
    end
end