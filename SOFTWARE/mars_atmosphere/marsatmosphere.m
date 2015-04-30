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
    end
    
    methods
        function obj = marsatmosphere()
            % General Mars Parameters
            % As given by the Mars fact sheet:
            % http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
            G = 6.67384e-11;                    %Nm2/kg2    Gravitational constant
            M_mars = 0.64174e24;                    %kg         Mars mass in kg
            V_mars = 16.318e19;                     %m^3        Mars Volume in m^3
            obj.r_mars = (V_mars*3/(4*pi))^(1/3);   %m          Mars Volumetric mean radius
            obj.GM_mars = G*M_mars;                 %Nm2/kg     G*M of mars for gravity
            
            % Atmosphere data from Mars-GRAM
            uid = fopen('input/atmosphere.txt');
            A = textscan(uid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
            fclose(uid);
            h_base = A{2}*1000;
            lat_base = A{3};
            lon_base = A{4};
            rho_base = A{5};
            T_base = A{6};
            CO2_base = A{13};
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
        
        function g = getg(obj, hq)
            g = obj.GM_mars ./ (hq+obj.r_mars).^2;
        end
    end
    
end