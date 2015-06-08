classdef aeroGeometry
% A definition of a geometry as used by modnewtonian, including helper functions

    properties
        %Geometry definition
        TriGeom;
        coords;
        tri;
        A_frontal;
        poly;
        
        %Calculated properties
        normals;
        areas;
        centers;
        centroid;
        combinations;
        directdistances;
        distancevectors;
        distances;
        
    end
    
    methods
        function obj = aeroGeometry(TriGeom, A_frontal, poly)
            obj.TriGeom = TriGeom;
            obj.A_frontal = A_frontal;
            obj.tri = TriGeom.ConnectivityList;
            obj.coords = TriGeom.Points';
            [obj.normals, obj.areas] = obj.calcNormalsAreas(obj.tri, obj.coords);
            obj.centers = obj.calcCellCenters(obj.tri, obj.coords);
            obj.centroid = obj.calcCentroid();
            obj.poly = poly;
%             [obj.distancevectors, obj.directdistances, obj.combinations, obj.distances] = obj.calcDirectDistances(obj.centers);

        end
        
        function [normals, areas] = calcNormalsAreas(~, tri, coords)
            %Calculate the surface normals of the given geometry.
            normals = zeros(size(tri'));
            areas = zeros(size(tri,1),1);
            for i = 1:size(tri,1)
                vec_a = coords(:,tri(i,2)) - coords(:,tri(i,1));
                vec_b = coords(:,tri(i,3)) - coords(:,tri(i,1));
                normals(:,i) = cross(vec_a, vec_b)/norm(cross(vec_a, vec_b));
                areas(i) = 0.5*norm(cross(vec_a, vec_b));
            end
        end
        
        function [] = plotPoly(obj)
            x = 0:0.01:1;
            plot(x, polyval(obj.poly, x)/sum(obj.poly));
        end
        
        function centers = calcCellCenters(~, tri, coords)
            % Calculate the cell centers of the given geometry
            centers = zeros(size(tri'));
            for i = 1:size(tri,1)
                vectors = [coords(:,tri(i,1)),coords(:,tri(i,2)), coords(:,tri(i,3))];
                centers(:,i) = mean(vectors,2);
            end
        end
        
        function centroid = calcCentroid(obj)
           centroid = obj.centers*obj.areas/sum(obj.areas);
        end
        
        function oppositePoint = getOpposite(obj, tbase, side)
            adjacents = obj.getAdjacentTriangles(tbase);
            temp = [3,1,2];
            for i = 1:3
                [baseside, checkside] = obj.getAdjacency(tbase, adjacents(i));
                if baseside == side
                    oppositePoint = obj.tri(adjacents(i),temp(checkside));
                    break;
                end
            end
        end
        
        function [baseside, checkside] = getAdjacency(obj, tbase, tcheck)
            % Base is the central triangle, check is the triangle of which
            % you want to know the connection
            combischeck = obj.tri(tcheck, :);
            combischeck = [combischeck(1),combischeck(2);combischeck(2),combischeck(3);combischeck(3),combischeck(1)];
            combisbase = obj.tri(tbase, :);
            combisbase = [combisbase(2),combisbase(1);combisbase(3),combisbase(2);combisbase(1),combisbase(3)];
            baseside = 0;
            checkside = 0;
            for basesidecount = 1:3
                for checksidecount = 1:3
                    if sum(combischeck(checksidecount,:)==combisbase(basesidecount,:))==2
                        baseside = basesidecount;
                        checkside = checksidecount;
                    end
                end
            end
        end
        
        function ind = getTriangle(obj, n1, n2, n3)
            % This function is not used as of now
            v = [n1; n2; n3];
            p = perms(v);
            ind = 0;
            for i = 1:size(p,1)
                equalrows = sum(obj.tri == ones(size(obj.tri)) * diag(p(i,:)), 2);
                if length(find(equalrows==3)) == 1
                    ind = find(equalrows==3);
                end
            end
        end
        
        function r = calcRadiusOfCurvature(obj, n1, n2, n3)
            A = obj.coords(:,n1);
            B = obj.coords(:,n2);
            C = obj.coords(:,n3);
            a = A-C;
            b = B-C;
            r = (norm(a)*norm(b)*norm(a-b))/(2*norm(cross(a,b)));
        end
        
        function indices = getAdjacentTriangles(obj, triangle)
            T = triangulation(obj.tri, obj.coords');
            indices = neighbors(T, triangle);
        end
        
        function points = getPointsOnXZPlane(obj, y)
            epsilon = 1e-10;
            points = find(obj.coords(2,:)>y-epsilon & obj.coords(2,:)<y+epsilon);
        end
        
        function points = getPointsOnYZPlane(obj, x)
            epsilon = 1e-10;            
            points = find(obj.coords(1,:)>x-epsilon & obj.coords(1,:)<x+epsilon);
        end
        
        function points = getPointsOnXYPlane(obj, z)
            epsilon = 1e-10;            
            points = find(obj.coords(3,:)>z-epsilon & obj.coords(3,:)<z+epsilon);
        end    
        
        function unconnectedTriangles = getUnconnectedTriangles(obj)
            unconnectedTriangles = [];
            for i = 1:size(obj.tri,1)
                if sum(isnan(obj.getAdjacentTriangles(i))) > 0
                    unconnectedTriangles(end+1) = i;
                end
            end
        end
        
        function [] = plotUnconnectedTriangles(obj, pauselength)
            unc = obj.getUnconnectedTriangles();
            for i = 1:length(unc)
                obj.plotTriangle(unc(i));
                pause(pauselength);
            end
        end
        
        function [] = plotPoint(obj, n)
        %PLOTPOINT Plot a point on the 3d plot
            plot3(obj.coords(1,n), obj.coords(2,n), obj.coords(3,n), 'o');
        end
        
        function [] = plotTriangle(obj, n)
        %PLOTTRIANGLE Plot the triangle
            t = obj.tri(n,:);
            obj.plotPoint(t(1));
            obj.plotPoint(t(2));
            obj.plotPoint(t(3));
        end     
   
        function triangles = getTrianglesOnPoint(obj, point)
            triangles = [find(obj.tri(:,1)==point);find(obj.tri(:,2)==point);find(obj.tri(:,3)==point)];
        end        
      
        
        function f = plotGeometry(obj, plotfaces, plotnormals)
            f = figure;
            h1 = axes;
            set(h1, 'Zdir', 'reverse');            
            hold on;
            caxis([0,1]);
            if plotfaces
                trisurf(obj.tri,obj.coords(1,:),obj.coords(2,:),obj.coords(3,:), rand(size(obj.tri(:,1))), 'EdgeColor', 'none');
            end
            axis equal;
            xlabel('x')
            ylabel('y')
            zlabel('z')
            if plotnormals
                quiver3(obj.centers(1,:), obj.centers(2,:), obj.centers(3,:), obj.normals(1,:), obj.normals(2,:), obj.normals(3,:))
            end
            zlim([min(obj.coords(3,:)), max(obj.coords(3,:))]);
            ylim([min(obj.coords(2,:)), max(obj.coords(2,:))]);
            view(40,30);
        end
        
        
        function f = plotValues(obj, values, variablename, colorrange, plotfaces, plotnormals)
            %f = plotValues(obj, values, variablename, colorrange, plotfaces, plotnormals)
            f = figure;
            h1 = axes;
            set(h1, 'Zdir', 'reverse');            
            hold on;
            if plotfaces
                trisurf(obj.tri,obj.coords(1,:),obj.coords(2,:),obj.coords(3,:), values, 'EdgeColor', 'none');
            end
            caxis(colorrange); %Cp goes from 0 to 2
            axis equal;
            h = colorbar;
            ylabel(h, variablename, 'interpreter', 'latex');
            xlabel('x')
            ylabel('y')
            zlabel('z')
            if plotnormals
                quiver3(obj.centers(1,:), obj.centers(2,:), obj.centers(3,:), obj.normals(1,:), obj.normals(2,:), obj.normals(3,:))
            end
            zlim([min(obj.coords(3,:)), max(obj.coords(3,:))]);
            ylim([min(obj.coords(2,:)), max(obj.coords(2,:))]);
            view(40,30);
            %Activate following code to plot vector of the wind
%             xlength = max(obj.coords(1,:))-min(obj.coords(1,:));
%             quiverV = - xlength * 0.5 * obj.V_array(:,end) / norm(obj.V_array(:,end));
%             quiverx = xlength*0.5 - quiverV(1) + max(obj.coords(1,:));
%             quiver3(quiverx,mean(obj.coords(2,:))-quiverV(2),mean(obj.coords(3,:))-quiverV(2),quiverV(1), quiverV(2), quiverV(3));
        end
        
        function [distances, obj] = getDistances(obj, triangle)
            distances = zeros(size(obj.tri,1),1);
            for i = 1:size(obj.tri,1)
                [distance, obj] = obj.getDistance(triangle,i);
                distances(i) = distance;
            end
        end

        function [distance, obj] = getDistance(obj, triangle1, triangle2)
            numcombinations = size(obj.combinations,1);
            index = find(sum(repmat(sort([triangle1, triangle2]),numcombinations,1)==obj.combinations,2)==2);
            if obj.distances(index) == 0
                distance = obj.calcDistance(triangle1, triangle2);
                obj.distances(index) = distance;
            elseif triangle1 == triangle2
                distance = 0;
            else
                distance = obj.distances(index);
            end
        end
        
        function [distancevectors, directdistances, combinations, distances] = calcDirectDistances(~, centers)
            combinations = nchoosek(1:size(centers,2),2);
            distancevectors = centers(:,combinations(:,2))-centers(:,combinations(:,1));
            directdistances = sqrt(sum(distancevectors.^2,1));
            distances = zeros(size(combinations,1),1);
        end
        
        function distance = calcDistance(obj, triangle1, triangle2)
            %calculate distances between points.
            numcombinations = size(obj.combinations,1);
            index = find(sum(repmat(sort([triangle1, triangle2]),numcombinations,1)==obj.combinations,2)==2);
            combination = obj.combinations(index,:);
            halfway = obj.centers(:,combination(1))+0.5*obj.distancevectors(:,index);
            pointdistances = sqrt(sum((obj.centers-repmat(halfway,1,size(obj.centers,2))).^2,1));
            [~,closestmidpoint] = min(pointdistances);
            if (closestmidpoint==combination(1))+(closestmidpoint==combination(2))>0
                distance = obj.directdistances(index);
            else
                distance1 = sum(obj.directdistances(sum(repmat(sort([combination(1),closestmidpoint]),numcombinations,1)==obj.combinations,2)==2));
                distance2 = sum(obj.directdistances(sum(repmat(sort([combination(2),closestmidpoint]),numcombinations,1)==obj.combinations,2)==2));
                distance = distance1 + distance2;
            end
        end
        
    end
    
    
    
end