
% General constants
G = 6.67384e-11;                    %Nm2/kg2    Gravitational constant

% As given by the Mars fact sheet:
% http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
M_mars = 0.64174e24;                %kg         Mars mass in kg
V_mars = 16.318e19;                 %m^3        Mars Volume in m^3
r_mars = (V_mars*3/(4*pi))^(1/3);   %m          Mars Volumetric mean radius
rho_mars = M_mars/V_mars;           %kg/m3      Mars Mean Density
g0_mars = 3.71;                     %m/s2       Mars surface acceleration
GM_mars = G*M_mars;                 %Nm2/kg     G*M of mars for gravity


