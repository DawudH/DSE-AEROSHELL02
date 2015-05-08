classdef modnewtonian
% Calculate aerodynamic properties of set geometry
    
    properties
        % Given at init
        a;
        gamma;
        tri;
        coords;
        center;
        
        % Given for calculation
        V_array;
        alpha_array;
        beta_array;
        M_array;
        
        % Calculated properties
        normals;
        cellcenters;
        areas;
        Cpmax_array;
        Cpdist_array;
        CRA_array;
        CMA_array;
    end
    
    methods
%         function obj = modnewtonian(TriGeom, gamma, a, center)
%             % Constructor: Set the geometry and gamma on Mars
%             obj.coords = TriGeom.Points';
%             obj.tri = TriGeom.ConnectivityList;
%             obj.gamma = gamma;
%             obj.a = a;
%             obj.center = center;
%             obj.normals = obj.calcsurfacenormals();
%             obj.cellcenters = obj.calccellcenters();
%         end
        
        function obj = modnewtonian(coords, tri, gamma, a, center)
            % Constructor: Set the geometry and gamma on Mars
            obj.coords = coords;
            obj.tri = tri;
            obj.gamma = gamma;
            obj.a = a;
            obj.center = center;
            [obj.normals, obj.areas] = obj.calcsurfacenormals();
            obj.cellcenters = obj.calccellcenters();
        end        
        
        function obj = calcAero(obj, V)
            % Calculate aerodynamic properties for the geometry for given V
            obj.M_array = [obj.M_array, sqrt(sum(V.^2))/obj.a];
            obj.Cpmax_array = [obj.Cpmax_array, obj.calcCp_max(obj.M_array(end), obj.gamma)];
            obj.V_array = [obj.V_array, V];
            obj.Cpdist_array = [obj.Cpdist_array, obj.calcCpdist()];
            obj.CRA_array = [obj.CRA_array, obj.calcForceCoeffs()];
        end
        
        function obj = calcAeroangle(obj, Vinf, alpha, beta)
            % Calculate aerodynamic properties with aoa and sideslip and V
            V = transformation(alpha, beta)*[Vinf;0;0];
            obj = obj.calcAero(V);
        end
        
        function CRA = calcForceCoeffs(obj)
            % Calculate aerodynamic force coefficients on the body
            CRA = - obj.normals * (obj.Cpdist_array(:,end) .* obj.areas);
        end
        
        function Cpdist = calcCpdist(obj)
            % Calculate CP-distribution based on modified newtonian
            Cpdist = zeros(size(obj.tri,1),1);
            Vhat = obj.V_array(:,end)/norm(obj.V_array(:,end));
            sinthetas = sum(obj.normals' * Vhat,2);
            sinthetas(sinthetas<0) = 0;
            Cpdist = obj.Cpmax_array(end)*sinthetas.^2;
        end
        
        
        function [SN, areas] = calcsurfacenormals(obj)
            %Calculate the surface normals of the given geometry.
            SN = zeros(size(obj.tri'));
            areas = zeros(size(obj.tri,1),1);
            for i = 1:size(obj.tri,1)
                vec_a = obj.coords(:,obj.tri(i,2)) - obj.coords(:,obj.tri(i,1));
                vec_b = obj.coords(:,obj.tri(i,3)) - obj.coords(:,obj.tri(i,1));
                SN(:,i) = cross(vec_a, vec_b)/norm(cross(vec_a, vec_b));
                areas(i) = norm(cross(vec_a, vec_b));
            end
        end
        
        function centers = calccellcenters(obj)
            % Calculate the cell centers of the given geometry
            centers = zeros(size(obj.tri'));
            for i = 1:size(obj.tri,1)
                vectors = [obj.coords(:,obj.tri(i,1)),obj.coords(:,obj.tri(i,2)), obj.coords(:,obj.tri(i,3))];
                centers(:,i) = mean(vectors,2);
            end
        end
        
        function Cp_max = calcCp_max(~, M, gamma)
            % Get the Cp_max as needed for modified newtonian theory
            Cp_max = 2./(gamma*M.^2).*((((gamma+1)^2*M.^2)/(4*gamma*M.^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M.^2)/(gamma+1))-1);
        end
        

    end
    
end

function T = transformation(alpha, beta)
    % Get the transformation matrix for angle of attack
    Ty = [[cos(-alpha), 0,          -sin(-alpha)];
          [0,           1,          0           ];
          [sin(-alpha), 0,          cos(-alpha) ]];

    Tz = [[cos(beta),   sin(beta),  0           ];
          [-sin(beta),  cos(beta),  0           ];
          [0,           0,          1           ]];
    T = Tz*Ty;
end

