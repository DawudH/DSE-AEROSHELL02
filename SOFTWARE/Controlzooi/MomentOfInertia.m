clear all
close all

%Moment coef
addpath('..\aerodynamic_coefficients');
aerocoef=aeroProperties();
aerocoef.getCML(2);

%Shape parameters
Di=4.5;
Do=12;
L=6;

Ds_cb=4;
m_cb=9000;

m_h=1000;
m_inf=150;

%Thruster parameters
Isp=234;
g0=9.81;
n=12

%Rotatian angles
mu_dd=5;
mu_d=20;
mu=60*ones(1,n);
alpha_dd=5;
alpha_d=20;
alpha_trim=10;
alpha_dive=0;
alpha=abs(alpha_trim-alpha_dive);

q=1000;
t_turn=10;

%Convert to radians
mu_dd=mu_dd/180*pi;
mu_d=mu_d/180*pi;
alpha_dd=alpha_dd/180*pi;
alpha_d=alpha_d/180*pi;
mu=mu/180*pi;
alpha=alpha/180*pi;

alpha_trim=alpha_trim/180*pi;
alpha_dive=alpha_dive/180*pi;


%Trim moment
CM_AL=abs(aerocoef.getCMLA(alpha_trim)-aerocoef.getCMLA(alpha_dive))
Myy_turn=CM_AL*q;


%------------------------------------
m_cbh=m_h-m_inf;
t_cb=Di-Ds_cb;

% %Compute MMOI
Ixx=1/4 * (0.5*m_cb*((Ds_cb+t_cb)^2+(Ds_cb-t_cb)^2) + 0.5*m_cbh*Di^2  +  0.5*m_inf+Do^2);
Iyy=1/12*(m_cb+m_h)*(3/4*((Ds_cb+t_cb)^2+(Ds_cb-t_cb)^2) +L^2);

%Rotational moment
Mxx=Ixx*mu_dd;
Myy=Iyy*alpha_dd;

%Compute thruster times
tx=2*min(sqrt(2*mu/mu_dd), mu_d/mu_dd);
ty=2*min(sqrt(2*alpha/alpha_dd), alpha_d/alpha_dd);

%Compute required thrust
Tx=Mxx/Di*2;
Ty=Myy/Di*2;
Ty_turn=Myy_turn/Di*2;

%Compute requried fuel
Massx=Tx*tx/Isp/g0;
Massx=sum(Massx)
Massy=Ty*ty/Isp/g0;
Massy_turn=Ty_turn*t_turn/Isp/g0;


% Delta V orbit
Isp=321;
t_burn=600;
delta_V=8.49;
delta_V_clean=3.49*3;
M0=10000;
g0=9.81;

%Massburn apocenter
M_p=M0-M0/(exp(delta_V/(Isp*g0)));
F=Isp*g0*M_p/t_burn
M_p_apo=M_p*2;

M_p_clean=M0-M0/(exp(delta_V_clean/(Isp*g0)));
F=Isp*g0*M_p_clean/t_burn


%Total mass without alpha control
M_tot=Massx+M_p_apo+M_p_clean;
V=M_tot/1.002*1.2;
M_tank=2.7086*10^-8 *V^3 -6.1703*10^-5 *V^2 +6.66290*10^-2 *V +1.3192;
M_thruster=8*1.6;

Mtank_fuel=M_tot+M_tank+M_thruster;

%Alternative
%Total mass with alpha control
M_tot_alpha=Massx+M_p_apo+Massy+Massy_turn+M_p_clean;

V=M_tot/1.002*1.2;
M_tank=2.7086*10^-8 *V^3 -6.1703*10^-5 *V^2 +6.66290*10^-2 *V +1.3192;
M_thruster=8*1.6;

Mtank_fuel_alpha=M_tot_alpha+M_tank+M_thruster;

table(M_p_apo,M_p_clean,M_tot,Mtank_fuel,Massy+Massy_turn,M_tot_alpha,Mtank_fuel_alpha)