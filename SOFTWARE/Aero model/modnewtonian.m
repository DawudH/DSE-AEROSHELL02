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
        Tmax_array;
        qmax_array;
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
            obj.CLCD_array = obj.CRA_aero_array(3,:)./obj.CRA_aero_array(1,:);
%             [T,q] = obj.calcHeatFlux();
%             obj.Tmax_array = [obj.Tmax_array, T];
%             obj.qmax_array = [obj.qmax_array, q];
            obj.CR_body_array = obj.CRA_body_array / obj.geom.A_frontal;
            obj.CR_aero_array = obj.CRA_aero_array / obj.geom.A_frontal;
            obj.CM_body_array = obj.CMA_body_array / obj.geom.A_frontal;
            obj.CM_aero_array = obj.CMA_aero_array / obj.geom.A_frontal;
        end
        
        function obj = calcAeroangle(obj, Vinf, alpha, beta)
            % Calculate aerodynamic properties with aoa and sideslip and V
            V = obj.Tba(alpha, beta)*[Vinf;0;0];
            obj.alpha_array = [obj.alpha_array, alpha];
            obj.beta_array = [obj.beta_array, beta];
            obj = obj.calcAero(V);
        end
        
        function CRAbody = calcForceCoeffsBody(obj)
            % Calculate aerodynamic force coefficients on the body
            CRAbody = - obj.geom.normals * (obj.Cpdist_array(:,end) .* obj.geom.areas);
        end
        
        function CRAaero = calcForceCoeffsAero(obj)
            CRAaero = diag([-1 1 -1])*(obj.Tab(obj.alpha_array(end), obj.beta_array(end))*obj.CRA_body_array(:,end));
        end
        
        function CMAbody = calcMomentCoeffsBody(obj)
            % Calculate aerodynamic moment coefficients on the body
            CMAbody = sum(cross(obj.geom.centers-repmat(obj.center',1,size(obj.geom.tri,1)), -obj.geom.normals * diag(obj.Cpdist_array(:,end) .* obj.geom.areas),1),2);
        end
        
        function CMAaero = calcMomentCoeffsAero(obj)
            CMAaero = obj.Tab(obj.alpha_array(end), obj.beta_array(end))*obj.CMA_body_array(:,end);
        end
        
        function Cpdist = calcCpdist(obj)
            % Calculate CP-distribution based on modified newtonian
            Vhat = obj.V_array(:,end)/norm(obj.V_array(:,end));
            sinthetas = sum(obj.geom.normals' * Vhat,2);
            sinthetas(sinthetas<0) = 0;
            Cpdist = obj.Cpmax_array(end)*sinthetas.^2;
        end

        function [Tmax, qmax] = calcHeatFlux(obj)
            % Get the stagnation heat flux and temperature
            M = 3;
            N = 0.5;
            Vinf = obj.a_inf*obj.M_array(end);
            Tmax = obj.T_inf*(obj.gamma-1)/2*obj.M_array(end)^2;
            
            [~, stagN] = max(obj.Cpdist_array(:,end));
            triangle = obj.geom.tri(stagN, :);
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
                radii = obj.geom.calcRadiusOfCurvature(point1, point2, 1);
            end
            qmax = 0;
            if max(radii) >=0
                qmax = obj.rho_inf^N*Vinf^M*1.83e-8*max(radii)^(-0.5);
            end
        end

        function obj = alphasweep(obj, Vinf, beta, alpha_start, alpha_end, dalpha)
            h = waitbar(0, 'Calculating...');
            for alpha = alpha_start:dalpha:alpha_end
                waitbar((alpha-alpha_start)/(alpha_end-alpha_start));
                obj = obj.calcAeroangle(Vinf, alpha, beta);
            end
            close(h);
        end
        
        function obj = betasweep(obj, Vinf, alpha, beta_start, beta_end, dbeta)
            for beta = beta_start:dbeta:beta_end
                obj = obj.calcAeroAngle(Vinf, alpha, beta);
            end
        end
        
        function obj = Msweep(obj, alpha, beta, M_start, M_end, dM)
            for M = M_start:M_end:dM
                obj = obj.calcAeroangle(obj.a_inf*M, alpha, beta);
            end
        end
        
        function obj = Vsweep(obj, alpha, beta, V_start, V_end, dV)
            for M = V_start:V_end:dV
                obj = obj.calcAeroangle(V, alpha, beta);
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

        function T = Tab(~, alpha, beta)
            % Get the transformation matrix for angle of attack
            Ty = [[cos(-alpha), 0,          -sin(-alpha)];
                  [0,           1,          0           ];
                  [sin(-alpha), 0,          cos(-alpha) ]];

            Tz = [[cos(beta),   sin(beta),  0           ];
                  [-sin(beta),  cos(beta),  0           ];
                  [0,           0,          1           ]];
            T = Tz*Ty;
        end

        function T = Tba(obj, alpha, beta)
            T = inv(obj.Tab(alpha, beta));
        end
        
        function obj = plotCp(obj, plotfaces, plotnormals)
            obj.geom.plotValues(obj.Cpdist_array(:,end), '$C_p$', [0, 2], plotfaces, plotnormals);
        end

        function obj = plots(obj, xarray, xlabeltxt, plotboolarray)
            % plotboolarray: ['cx', 'cy', 'cz', 'cd', 'cl', 'cs']
            for i = plotboolarray
                j = cell2mat(i{1});
                switch j
                    case 'cla'
                        figure;
                        plot(xarray, obj.CRA_aero_array(3,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'csa'
                        figure;
                        plot(xarray, obj.CRA_aero_array(2,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cda'
                        figure;
                        plot(xarray, obj.CRA_aero_array(1,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cxa'
                        figure;
                        plot(xarray, obj.CRA_body_array(1,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cya'
                        figure;
                        plot(xarray, obj.CRA_body_array(2,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cza'
                        figure;
                        plot(xarray, obj.CRA_body_array(3,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'clcd'
                        figure;
                        plot(xarray, obj.CLCD_array);
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cmya'
                        figure;
                        plot(xarray, obj.CMA_aero_array(2,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cmycx'
                        figure;
                        plot(xarray, obj.CMA_body_array(2,:)./obj.CRA_body_array(1,:));
                        ylabel(j);
                        xlabel(xlabeltxt);                     
                    case 'q'
                        figure;
                        plot(xarray, obj.qmax_array);
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'T'
                        figure;
                        plot(xarray, obj.Tmax_array);
                        ylabel(j);
                        xlabel(xlabeltxt);
                    otherwise
                        disp(strcat('Not a valid input for plotting: ', num2str(i)));
                end
            end
        end
    end
end

