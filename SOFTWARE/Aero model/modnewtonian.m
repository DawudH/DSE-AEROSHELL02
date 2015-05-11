classdef modnewtonian
% Calculate aerodynamic properties of set geometry
    
    properties
        % Given at init
        a_inf;
        gamma;
        tri;
        coords;
        center;
        rho_inf;
        T_inf;
        
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
        CRA_body_array;
        CRA_aero_array;
        CMA_body_array;
        CMA_aero_array;
        CLCD_array;
        Tmax_array;
        qmax_array;
    end
    
    methods
%         function obj = modnewtonian(TriGeom, gamma, a, center)
%             % Constructor: Set the geometry and gamma on Mars
%             obj.coords = TriGeom.Points';
%             obj.tri = TriGeom.ConnectivityList;
%             obj.gamma = gamma;
%             obj.a_inf = a;
%             obj.center = center;
%             obj.normals = obj.calcsurfacenormals();
%             obj.cellcenters = obj.calccellcenters();
%         end
        
        function obj = modnewtonian(coords, tri, gamma, a, center, rho, T)
            % Constructor: Set the geometry and gamma on Mars
            obj.coords = coords;
            obj.tri = tri;
            obj.gamma = gamma;
            obj.a_inf = a;
            obj.center = center;
            [obj.normals, obj.areas] = obj.calcsurfacenormals();
            obj.cellcenters = obj.calccellcenters();
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
            [T,q] = obj.calcHeatFlux();
            obj.Tmax_array = [obj.Tmax_array, T];
            obj.qmax_array = [obj.qmax_array, q];
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
            CRAbody = - obj.normals * (obj.Cpdist_array(:,end) .* obj.areas);
        end
        
        function CRAaero = calcForceCoeffsAero(obj)
            CRAaero = diag([-1 1 -1])*(obj.Tab(obj.alpha_array(end), obj.beta_array(end))*obj.CRA_body_array(:,end));
        end
        
        function CMAbody = calcMomentCoeffsBody(obj)
            % Calculate aerodynamic moment coefficients on the body
            CMAbody = sum(cross(obj.cellcenters, -obj.normals * diag(obj.Cpdist_array(:,end) .* obj.areas),1),2);
        end
        
        function CMAaero = calcMomentCoeffsAero(obj)
            CMAaero = obj.Tab(obj.alpha_array(end), obj.beta_array(end))*obj.CMA_body_array(:,end);
        end
        
        function Cpdist = calcCpdist(obj)
            % Calculate CP-distribution based on modified newtonian
            Vhat = obj.V_array(:,end)/norm(obj.V_array(:,end));
            sinthetas = sum(obj.normals' * Vhat,2);
            sinthetas(sinthetas<0) = 0;
            Cpdist = obj.Cpmax_array(end)*sinthetas.^2;
        end
        
        function r = radiusOfCurvature(obj, n1, n2, n3)
            A = obj.coords(:,n1);
            B = obj.coords(:,n2);
            C = obj.coords(:,n3);
            a = A-C;
            b = B-C;
            r = (norm(a)*norm(b)*norm(a-b))/(2*norm(cross(a,b)));
        end
        
        function [Tmax, qmax] = calcHeatFlux(obj)
            % Get the stagnation heat flux and temperature
            M = 3;
            N = 0.5;
            Vinf = obj.a_inf*obj.M_array(end);
            Tmax = obj.T_inf*(1+obj.gamma)/2*obj.M_array(end)^2;
            
            [cp, stagN] = max(obj.Cpdist_array(:,end));
            maxtri = max(max(obj.tri));
            
            points = perms(obj.tri(stagN,:));
            combis = points(:,1:2);
            nextpoints = 2*combis(:,2)-combis(:,1);
            nextpoints = nextpoints + maxtri;
            nextpoints = mod(nextpoints, maxtri);
            radii = zeros(size(nextpoints));
            for i = 1:length(nextpoints)
                radius = obj.radiusOfCurvature(combis(i,1), combis(i,2), nextpoints(i));
                radii(i) = radius;
            end
            qmax = obj.rho_inf^0.5*Vinf^3*1.83e-8*max(radii)^(-0.5);
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

        function [SN, areas] = calcsurfacenormals(obj)
            %Calculate the surface normals of the given geometry.
            SN = zeros(size(obj.tri'));
            areas = zeros(size(obj.tri,1),1);
            for i = 1:size(obj.tri,1)
                vec_a = obj.coords(:,obj.tri(i,2)) - obj.coords(:,obj.tri(i,1));
                vec_b = obj.coords(:,obj.tri(i,3)) - obj.coords(:,obj.tri(i,1));
                SN(:,i) = cross(vec_a, vec_b)/norm(cross(vec_a, vec_b));
                areas(i) = 0.5*norm(cross(vec_a, vec_b));
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
        
        function obj = plotCp(obj, plotfaces, plotnormals)
            figure;
            h1 = axes;
            set(h1, 'Zdir', 'reverse');            
            hold on;
            if plotfaces
                trisurf(obj.tri,obj.coords(1,:),obj.coords(2,:),obj.coords(3,:), obj.Cpdist_array(:,end));
            end
            caxis([0,2]); %Cp goes from 0 to 2
            axis equal;
            colorbar;
            xlabel('x')
            ylabel('y')
            zlabel('z')
            if plotnormals
                quiver3(obj.cellcenters(1,:), obj.cellcenters(2,:), obj.cellcenters(3,:), obj.normals(1,:), obj.normals(2,:), obj.normals(3,:))
            end
            xlength = max(obj.coords(1,:))-min(obj.coords(1,:));
            quiverV = - xlength * 0.5 * obj.V_array(:,end) / norm(obj.V_array(:,end));
            quiverx = xlength*0.5 - quiverV(1);
            quiver3(quiverx,mean(obj.coords(2,:))-quiverV(2),mean(obj.coords(3,:))-quiverV(2),quiverV(1), quiverV(2), quiverV(3));

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
        
        function obj = plots(obj, xarray, xlabeltxt, plotboolarray)
            % plotboolarray: ['cx', 'cy', 'cz', 'cd', 'cl', 'cs']
            for i = plotboolarray
                j = cell2mat(i{1});
                switch j
                    case 'cl'
                        figure;
                        plot(xarray, obj.CRA_aero_array(3,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cs'
                        figure;
                        plot(xarray, obj.CRA_aero_array(2,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cd'
                        figure;
                        plot(xarray, obj.CRA_aero_array(1,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cx'
                        figure;
                        plot(xarray, obj.CRA_body_array(1,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cy'
                        figure;
                        plot(xarray, obj.CRA_body_array(2,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cz'
                        figure;
                        plot(xarray, obj.CRA_body_array(3,:));
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'clcd'
                        figure;
                        plot(xarray, obj.CLCD_array);
                        ylabel(j);
                        xlabel(xlabeltxt);
                    case 'cmy'
                        figure;
                        plot(xarray, obj.CMA_aero_array(2,:));
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

