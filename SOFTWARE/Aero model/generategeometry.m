function [ coords, tri, A ] = generategeometry( type, q )
%[ coords, tri, A ] = generategeometry( type, q )
%   Generate used geometries


    switch type
        case 'horizontalplate'
            xvector = [0 1 0 1];
            yvector = [0 0 1 1];
            zvector = [0 0 0 0];
            coords = [xvector;yvector;zvector];
            tri = [1 2 3;2 4 3];
            A = 1;
        case 'verticalplate'
            xvector = [0 0 0 0];
            yvector = [0 0 1 1];
            zvector = [0 1 1 0];
            coords = [xvector;yvector;zvector];        
            tri = [1 4 2;2 4 3];
            A = 1;
        case 'sphere12m'
            TriGeom = TriMeshGen(q, 6, 6, 6, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList;
            A = pi*6^2;
        case 'pastille12m1.5m'
            TriGeom = TriMeshGen(q, 6, 1.5, 6, 's');
            coords = TriGeom.Points';
            tri = TriGeom.ConnectivityList; 
            A = pi*6^2;
        otherwise
            warning(strcat('The following type is not supported: ', num2str(type)));
            xvector = [0,0,0,0,-1,-1];
            yvector = [0,1,1,0,1,0];
            zvector = [0,0,1,1,0,0];
            coords = [xvector;yvector;zvector;];
            tri = [1 2 4;2 3 4;1 5 2;1 6 5];
            A = 1;
    end
end