% Parameter file as general start - contains general constants and Mars
% parameters

% General constants
G = 6.67384e-11;                    %Nm2/kg2    Gravitational constant


% General Mars Parameters
% As given by the Mars fact sheet:
% http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
% Most values are calculated, agree with the values given by the fact sheet
M_mars = 0.64174e24;                %kg         Mars mass in kg
V_mars = 16.318e19;                 %m^3        Mars Volume in m^3
r_mars = (V_mars*3/(4*pi))^(1/3);   %m          Mars Volumetric mean radius
rho_mars = M_mars/V_mars;           %kg/m3      Mars Mean Density
