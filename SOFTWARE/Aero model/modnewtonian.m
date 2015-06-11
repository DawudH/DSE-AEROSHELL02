classdef modnewtonian
% Calculate aerodynamic properties of set geometry
    
    properties
        % Given at init
        a_inf;
        gamma;
        center;
        rho_inf;
        T_inf;
        geom
        
        % Given for calculation
        V_array;
        alpha_array;
        beta_array;
        phi_array;
        M_array;
        
        % Calculated properties
        Cpmax_array;
        Cpdist_array;
        CRA_body_array;
        CRA_aero_array;
        CMA_body_array;
        CMA_aero_array;
        CR_body_array;
        CR_aero_array;
        CM_body_array;
        CM_aero_array;
        CLCD_array;
        CG_offset;
        Cmyalpha
    end
    
    methods
        
        function obj = modnewtonian(geom, gamma, a, center, rho, T)
            % Constructor: Set the geometry and gamma on Mars
            obj.geom = geom;
            obj.gamma = gamma;
            obj.a_inf = a;
            obj.center = center;
            obj.rho_inf = rho;
            obj.T_inf = T;
        end
        
        function obj = calcAero(obj, V)
            % Calculate aerodynamic properties for the geometry for given V
            obj.M_array = [obj.M_array, sqrt(sum(V.^2))/obj.a_inf];
            obj.Cpmax_array = [obj.Cpmax_array, obj.calcCp_max(obj.M_array(end), obj.gamma)];
            obj.V_array = [obj.V_array, V];
            obj.Cpdist_array = [obj.Cpdist_array, obj.calcCpdist()];
            obj.CRA_body_array = [obj.CRA_body_array, obj.calcForceCoeffsBody()];
            obj.CRA_aero_array = [obj.CRA_aero_array, obj.calcForceCoeffsAero()];
            obj.CMA_body_array = [obj.CMA_body_array, obj.calcMomentCoeffsBody()];
            obj.CMA_aero_array = [obj.CMA_aero_array, obj.calcMomentCoeffsAero()];
            obj = obj.calcPostCalculations();
        end
        
        function obj = calcPostCalculations(obj)
            obj.CLCD_array = obj.CRA_aero_array(3,:)./obj.CRA_aero_array(1,:);
            obj.CR_body_array = obj.CRA_body_array / obj.geom.A_frontal;
            obj.CR_aero_array = obj.CRA_aero_array / obj.geom.A_frontal;
            obj.CM_body_array = obj.CMA_body_array / obj.geom.A_frontal;
            obj.CM_aero_array = obj.CMA_aero_array / obj.geom.A_frontal;
            obj.CG_offset = obj.CM_body_array(2,:)./obj.CR_body_array(1,:);
            obj.Cmyalpha = obj.calcCmyalpha();
        end
        
        function obj = calcAeroangle(obj, Vinf, alpha, beta, phi)
            % Calculate aerodynamic properties with aoa and sideslip and V
            V = obj.Tba(alpha, beta, phi)*[Vinf;0;0];
            obj.alpha_array = [obj.alpha_array, alpha];
            obj.beta_array = [obj.beta_array, beta];
            obj.phi_array = [obj.phi_array, phi];
            obj = obj.calcAero(V);
        end
        
        function CRAbody = calcForceCoeffsBody(obj)
            % Calculate aerodynamic force coefficients on the body
            CRAbody = - obj.geom.normals * (obj.Cpdist_array(:,end) .* obj.geom.areas);
        end
        
        function CRAaero = calcForceCoeffsAero(obj)
            CRAaero = diag([-1 1 -1])*(obj.Tab(obj.alpha_array(end), obj.beta_array(end), obj.phi_array(end))*obj.CRA_body_array(:,end));
        end
        
        function CMAbody = calcMomentCoeffsBody(obj)
            % Calculate aerodynamic moment coefficients on the body
            CMAbody = sum(cross(obj.geom.centers-repmat(obj.center',1,size(obj.geom.tri,1)), -obj.geom.normals * diag(obj.Cpdist_array(:,end) .* obj.geom.areas),1),2);
        end
        
        function CMAaero = calcMomentCoeffsAero(obj)
            CMAaero = obj.Tab(obj.alpha_array(end), obj.beta_array(end), obj.phi_array(end))*obj.CMA_body_array(:,end);
        end
        
        function Cpdist = calcCpdist(obj)
            % Calculate CP-distribution based on modified newtonian
            Vhat = obj.V_array(:,end)/norm(obj.V_array(:,end));
            sinthetas = sum(obj.geom.normals' * Vhat,2);
            sinthetas(sinthetas<0) = 0;
            Cpdist = obj.Cpmax_array(end)*sinthetas.^2;
        end
        
        function CmyAalpha = calcCmyAalpha(obj)
            CmyAalpha = diff(obj.CMA_body_array(2,:))./diff(obj.alpha_array);
        end
        
        function Cmyalpha = calcCmyalpha(obj)
            Cmyalpha = diff(obj.CM_body_array(2,:))./diff(obj.alpha_array);
        end
        
        
        
        function [Tmax_array, qmax_array] = calcStagnationHeatFluxes(obj)
            Tmax_array = zeros(size(obj.alpha_array));
            qmax_array = zeros(size(obj.alpha_array));
                        
            % Get the stagnation heat flux and temperature
            for j = 1:length(obj.alpha_array)
%             Vhat = obj.V_array(:,end)/norm(obj.V_array(:,end));
                Vinf = norm(obj.V_array(:,j));
    %             sinthetas = obj.geom.normals' * Vhat;
    %             sinthetas(sinthetas<0) = 0;
    %             costhetas = cos(asin(sinthetas));

                Tmax = obj.T_inf*(obj.gamma-1)/2*obj.M_array(j)^2;
                [~, stagN] = max(obj.Cpdist_array(:,j));
                triangle = obj.geom.tri(stagN, :);

                %calculate radius of curvature
                if sum(triangle==[1 1 1])==0
                    opposites = zeros(1,3);
                    for i = 1:3
                        opposites(i) = obj.geom.getOpposite(stagN, i);
                    end

                    checkMatrix = [opposites(3),triangle(1), opposites(1); ...
                                    opposites(1), triangle(2), opposites(2); ...
                                    opposites(2), triangle(3), opposites(3)];

                    radii = zeros(3,1);
                    for i = 1:length(radii)
                        radii(i) = obj.geom.calcRadiusOfCurvature(checkMatrix(i,1), checkMatrix(i,2), checkMatrix(i,3));
                    end
                    radius = max(radii);
                else %If point on centerpoint fuck
                    trianglesincircle = sum(sum(obj.geom.tri == ones(size(obj.geom.tri))));
                    point1 = triangle(1);
                    if point1 == 1
                        point1 = triangle(2);
                    end
                    overflowvector = [0.5*trianglesincircle+1:trianglesincircle, 1:0.5*trianglesincircle];
                    point2 = overflowvector(point1);
                    if obj.geom.coords(2, point1) - obj.geom.coords(2,point2) < 1e-12
                        point2 = point2 + 1;
                    end
                    radius = obj.geom.calcRadiusOfCurvature(point1, point2, 1);
                end

                %stagnation point values

                qmax = 0;
                if radius >=0
                    M = 3;
                    N = 0.5;
                    C = 1.83e-8*max(radius)^-.5;                
                    qmax = obj.rho_inf ^ N * Vinf ^ M * C;
                end
                qmax_array(j) = qmax;
                Tmax_array(j) = Tmax;
            end

    %             %distribution values
    %             xt = obj.geom.getDistances(stagN);
    %             
    %             N = 0.8;
    %             if Vinf <= 3962
    %                 M = 3.37;
    %                 C = (3.89e-8)*(costhetas.^1.78).*(sinthetas.^1.6).*(xt.^-.2).*(Tw/556).^(-.25).*(1-1.11*Tw/Tmax);
    %             else
    %                 M = 3.7;
    %                 C = (2.2e-9)*(costhetas.^2.08).*(sinthetas.^1.6).*(xt.^-.2).*(1-1.11*Tw/Tmax);
    %             end
    %             qw = obj.rho_inf(end) ^ N * Vinf ^ M * C;
    %             qw = sinthetas*qmax;
        end

        function obj = alphasweep(obj, Vinf, beta, phi, alpha_start, alpha_end, dalpha)
%             h = waitbar(0, 'Calculating...');
            for alpha = alpha_start:dalpha:alpha_end
%                 waitbar((alpha-alpha_start)/(alpha_end-alpha_start));
                V = obj.Tba(alpha, beta, phi)*[Vinf;0;0];
                obj.alpha_array = [obj.alpha_array, alpha];
                obj.beta_array = [obj.beta_array, beta];
                obj.phi_array = [obj.phi_array, phi];
                obj = obj.calcAero(V);
            end
%             close(h);
            obj.calcPostCalculations();
        end
        
        function obj = betasweep(obj, Vinf, alpha, phi, beta_start, beta_end, dbeta)
            for beta = beta_start:dbeta:beta_end
                obj = obj.calcAeroAngle(Vinf, alpha, beta, phi);
            end
        end
        
        function obj = phisweep(obj, Vinf, alpha, beta, phi_start, phi_end, dphi)
            for phi = phi_start:dphi:phi_end
                obj = obj.calcAeroAngle(Vinf, alpha, beta, phi);
            end
        end
        
        function obj = Msweep(obj, alpha, beta, phi, M_start, M_end, dM)
            for M = M_start:M_end:dM
                obj = obj.calcAeroangle(obj.a_inf*M, alpha, beta, phi);
            end
        end
        
        function obj = Vsweep(obj, alpha, beta, phi, V_start, V_end, dV)
            for M = V_start:V_end:dV
                obj = obj.calcAeroangle(V, alpha, beta, phi);
            end
        end
        
        function Cp_max = calcCp_max(~, M, gamma)
            % Get the Cp_max as needed for modified newtonian theory
            Cp_max = 2./(gamma*M.^2).*((((gamma+1)^2*M.^2)/(4*gamma*M.^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M.^2)/(gamma+1))-1);
        end

        function cp = calcCpOnPoint(obj, point)
            triangles = obj.getTrianglesOnPoint(point);
            cp = mean(obj.Cpdist_array(triangles, end));
        end

        function T = Tab(~, alpha, beta, phi)
            % Get the transformation matrix for angle of attack
            Tx = [[1,           0,          0           ];
                  [0,           cos(phi),   sin(phi)    ];
                  [0,           -sin(phi),  cos(phi)    ]];
              
            Ty = [[cos(-alpha), 0,          -sin(-alpha)];
                  [0,           1,          0           ];
                  [sin(-alpha), 0,          cos(-alpha) ]];

            Tz = [[cos(beta),   sin(beta),  0           ];
                  [-sin(beta),  cos(beta),  0           ];
                  [0,           0,          1           ]];
            T = Tz*Ty*Tx;
        end

        function T = Tba(obj, alpha, beta, phi)
            T = inv(obj.Tab(alpha, beta, phi));
        end
        
        function obj = plotCp(obj, plotfaces, plotnormals)
            obj.geom.plotValues(obj.Cpdist_array(:,end), '$C_p$', [0, 2], plotfaces, plotnormals);
        end

        function obj = plots(obj)
            obj.geom.plotGeometry(true, false);
            figure;
            plot(rad2deg(obj.alpha_array), obj.CR_aero_array');
            xlabel('alpha')
            ylabel('CR');
            legend('CD', 'CS', 'CL');
            
            figure;
            plot(rad2deg(obj.alpha_array), obj.CR_body_array');
            xlabel('alpha')
            ylabel('CR');
            legend('CX', 'CY', 'CZ');
            
            figure;
            plot(rad2deg(obj.alpha_array), obj.CM_aero_array');
            xlabel('alpha')
            ylabel('CM');
            legend('CMx', 'CMy', 'CMZ');
            
            figure; %L over D
            plot(rad2deg(obj.alpha_array), obj.CLCD_array);
            xlabel('alpha');
            ylabel('CL/CD');
            legend('CL/CD');
            
            figure;
            plot(rad2deg(obj.alpha_array), obj.CM_body_array(2,:)./obj.CR_body_array(1,:));
            xlabel('alpha');
            ylabel('CG offset [m]');
            legend('CG offset [m]');
            
            figure;
            plot(obj.CLCD_array, obj.CM_body_array(2,:)./obj.CR_body_array(1,:));
            xlabel('CL/CD');
            ylabel('CG offset [m]');
            legend('CG offset [m]');
            
        end
    end
end

