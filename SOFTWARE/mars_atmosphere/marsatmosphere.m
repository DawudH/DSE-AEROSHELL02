classdef marsatmosphere
%atmospheric calculations of Mars. height in m
    
    properties
        latmesh;
        lonmesh;
        hmesh;
        Tmesh;
        rhomesh;
        GM_mars;
        hmax;
        r_mars;
        CO2_base;
        gamma_mars;
        Rm_mars;
        cheaprho;
        cheapT;
        cheaph;
        cheapLat;
        cheapLon;
        Tinterpolant;
        rhointerpolant;
    end
    
    methods
        function obj = marsatmosphere()
            % General Mars Parameters
            % As given by the Mars fact sheet:
            % http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
            G = 6.67384e-11;                    %Nm2/kg2    Gravitational constant
            M_mars = 0.64174e24;                    %kg         Mars mass in kg
            obj.r_mars = 3396000;   %m          Mars Volumetric mean radius
            obj.GM_mars = G*M_mars;                 %Nm2/kg     G*M of mars for gravity
            
            % Atmosphere data from Mars-GRAM
            location = which('atmosphere.txt');
            uid = fopen(location);
            A = textscan(uid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
            fclose(uid);
            h_base = A{2}*1000;
            lat_base = A{3};
            lon_base = A{4};
            rho_base = A{5};
            T_base = A{6};
            obj.CO2_base = A{13};
            N2_base = A{14};
            Ar_base = A{15};
            O2_base = A{16};
            O_base = A{17};
            He_base = A{18};
            H2_base = A{19};
            H_base = A{20};
            H2O_base = A{21};
            
            obj.hmax = max(h_base);
            
            latnom = length(unique(lat_base));
            lonnom = length(unique(lon_base));
            hnom = length(unique(h_base));

            obj.Tmesh = reshape(T_base,[lonnom, latnom, hnom]);
            obj.rhomesh = reshape(rho_base,[lonnom, latnom, hnom]);

            [obj.latmesh,obj.lonmesh,obj.hmesh] = meshgrid(unique(lat_base), unique(lon_base), unique(h_base));
            
            obj.gamma_mars = 1.29;
            obj.Rm_mars = 191.8;
            
            obj.cheapLat = 0;
            obj.cheapLon = 180;
            
            obj.cheaph = 0:5000:4e5;
            obj.cheaprho = obj.getDensity(obj.cheapLat*zeros(size(obj.cheaph)), obj.cheapLon*zeros(size(obj.cheaph)), obj.cheaph);
            obj.cheapT = obj.getTemperature(obj.cheapLat*zeros(size(obj.cheaph)), obj.cheapLon*zeros(size(obj.cheaph)), obj.cheaph);
            obj.rhointerpolant = griddedInterpolant(obj.cheaph, obj.cheaprho);
            obj.Tinterpolant = griddedInterpolant(obj.cheaph, obj.cheapT);
        end
        
        function rho = getCheapDensity(obj, hq)
            if hq > obj.hmax
                rho = 0;
            else
                rho = obj.rhointerpolant(hq);
%                 rho = interp1(obj.cheaph, obj.cheaprho, hq, 'linear');
            end            
        end
        
        function T = getCheapTemperature(obj, hq)
            if hq > obj.hmax
                T = obj.cheapT(end);
            else
                T = obj.Tinterpolant(hq);
%                 T = interp1(obj.cheaph, obj.cheapT, hq, 'nearest');
            end            
        end        
        
        function rho = getDensity(obj, latq, lonq, hq) %only for scalars
            if latq < -90 % if latitude lower than 90 degrees
                latq  = -180-latq;
            end
            if latq > 90
                latq = 180-latq;
            end
            lonq = mod(lonq,360);
            if hq > obj.hmax
                rho = 0;
            else
                rho = interp3(obj.latmesh, obj.lonmesh, obj.hmesh, obj.rhomesh, latq,lonq,hq,'cubic');
            end
        end
        
        function T = getTemperature(obj, latq, lonq, hq) %only for scalars
            if latq < -90 % if latitude lower than 90 degrees
                latq  = -180-latq;
            end
            if latq > 90
                latq = 180-latq;
            end
            lonq = mod(lonq,360);
            if hq > obj.hmax
                hq = obj.hmax;
            end
            T = interp3(obj.latmesh, obj.lonmesh, obj.hmesh, obj.Tmesh, latq,lonq,hq,'cubic');
        end
        
        function a = getSpeedofsound(obj, latq, lonq, hq)
            T = obj.getTemperature(latq, lonq, hq);
            a = sqrt(obj.gamma_mars*obj.Rm_mars*T);
        end
        
        function g = getg(obj, hq)
            g = obj.GM_mars ./ (hq+obj.r_mars).^2;
        end
    end
    
end