x_1=0.3;
x_2=1.3;
z_cg=linspace(1.165,4.,100);

phi_1=atan(x_1./z_cg);
phi_2=atan(x_2./z_cg);

plot(z_cg,rad2deg(phi_2-phi_1));
xlabel('Z_{cg}');
ylabel('\Delta \phi');




theta=45; % surface inclination in degrees, w.r.t. body-symmetric axis (body frame)
alpha=0; % angle of attack
              
q=2000; % Dynamic pressure
gamma=1.4;
M=30;
num_flaps=2;
M_trim_deg=1717;
flap_arm=6;
flap_area=1; % m^2, per flap
flap_centroid=[4, 0, 6] % w.r.t. C.G.

Cp_max = 2/(gamma*M^2)*((((gamma+1)^2*M^2)/(4*gamma*M^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M^2)/(gamma+1))-1);
flap_normal=[-cos(deg2rad(theta+alpha)),0,-sin(deg2rad(theta+alpha))]; % in body frame
CR_A=-Cp_max*(sin(deg2rad(90+theta+alpha))^2).*flap_area*flap_normal;

CM_A_flap=cross(flap_centroid,CR_A)