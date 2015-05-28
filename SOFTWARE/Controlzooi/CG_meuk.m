x_1=0.3;
x_2=1.3;
z_cg=linspace(1.165,4.,100);

phi_1=atan(x_1./z_cg);
phi_2=atan(x_2./z_cg);

plot(z_cg,rad2deg(phi_2-phi_1));
xlabel('Z_{cg}');
ylabel('\Delta \phi');

theta=deg2rad(90);

q=2000; % Dynamic pressure
gamma=1.4;
M=30;
Cp_max = 2./(gamma*M^2).*((((gamma+1)^2*M^2)/(4*gamma*M^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*M^2)/(gamma+1))-1);
num_flaps=2;
M_trim_deg=1717;
flap_arm=6;
Flap_area=1; % m^2, per flap

flap_normal=[-cos(theta),0,sin(theta)];
[C_x,C_y,C_z]=-Cp_max*(sin(theta)^2)*Flap_area

F_trim_deg=M_trim_deg/flap_arm;
F_trim_deg_perflap=F_trim_deg/num_flaps;
CD_trim_deg_flap=F_trim_deg_perflap/(q*Flap_area);
