classdef modnewtonian
% Calculate aerodynamic properties of set geometry
    
    properties
        % Given at init
        a;
        gamma;
        geometry;
        center;
        
        % Given for calculation
        V_array;
        alpha_array;
        beta_array;
        M_array;
        
        % Calculated properties
        Cpdistr_array;
        CLS_array;
        CDS_array;
        CMS_array;
    end
    
    methods
        function obj = modnewtonian(geometry, gamma, a, center)
            % Constructor: Set the geometry and gamma on Mars
            obj.geometry = geometry;
            obj.gamma = gamma;
            obj.a = a;
            obj.center = center;
        end
        
        function [CLS, CDS, CMS] = calcAero(obj, V)
            % Calculate aerodynamic properties for the geometry for given V
            obj.M = sqrt(sum(V.^2))/obj.a;
            obj.Cpdist = [obj.Cpdist; getCpdist(obj)];
            
            % Calculate!
            
        end
        
        function [CLS, CDS, CMS] = calcAeroangle(obj, Vinf, alpha, beta)
            % Calculate aerodynamic properties with aoa and sideslip and V
            V = transformation(alpha, beta)*[0,0,Vinf];
            [CLS, CDS, CMS] = calcAero(obj, V);
        end
        
        function Cpdist = calcCpdist(obj)
            % Calculate CP-distribution based on modified newtonian
            Cpdist = [];            
        end
        
        function CL = calcCL
        
        function Cp_max = getCp_max(obj, M)
            % Get the Cp_max as needed for modified newtonian theory
            gamma = obj.gamma; %To make the calculation shorter
            Cp_max = 2./(gamma*M.^2).*((((gamma+1)^2*M.^2)/(4*gamma*M.^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M.^2)/(gamma+1))-1);
            obj.Cp_max = Cp_max;
        end
        
        function T = transformation(alpha, beta)
            % Get the transformation matrix for angle of attack
            Ty = [[cos(-alpha), 0,  -sin(-alpha)];
                  [0,           1,  0           ];
                  [sin(-alpha), 0,  cos(-alpha) ]];
            Tz = [[cos(beta),   sin(beta),  0];
                  [-sin(beta),  cos(beta),  0];
                  [0,           0,          1]];
            T = Tz*Ty;
        end
        
        
        
        
        
    end
    
end