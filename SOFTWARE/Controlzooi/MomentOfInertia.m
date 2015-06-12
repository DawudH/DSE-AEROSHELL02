clear all
close all

addpath('..\aerodynamic_coefficients');
aerocoef=aeroProperties();
aerocoef.getCML(2);

Di=4.5;
Do=12;
L=6;

Ds_cb=3.5;
m_cb=9000;

m_h=1000;
m_inf=150;

Isp=234;
g0=9.81;

mu_dd=5;
mu=41;
alpha_dd=5;
alpha_trim=5;
alpha_dive=0;
alpha=abs(alpha_trim-alpha_dive)

mu_dd=mu_dd/180*pi;
alpha_dd=alpha_dd/180*pi;
mu=mu/180*pi;
alpha=alpha/180*pi;

alpha_trim=alpha_trim/180*pi;
alpha_dive=alpha_dive/180*pi;

q=1000;

CM_AL=abs(aerocoef.getCMLA(alpha_trim)-aerocoef.getCMLA(alpha_dive))
Myy_turn=CM_AL*q;
t_turn=5;

m_cbh=m_h-m_inf;
t_cb=Di-Ds_cb;

Ixx=1/4 * (0.5*m_cb*((Ds_cb+t_cb)^2+(Ds_cb+t_cb)^2) + 0.5*m_cbh*Di^2  +  0.5*m_inf+Do^2);
Iyy=1/12*(m_cb+m_h)*(3/4*((Ds_cb+t_cb)^2+(Ds_cb+t_cb)^2) +L^2);


Mxx=Ixx*mu_dd*2;
Myy=Iyy*alpha_dd*2;


txx=2*sqrt(2*mu/mu_dd);
tyy=2*sqrt(2*alpha/alpha_dd);

Txx=Mxx/Di;
Tyy=Myy/Di;
Tyy_turn=Myy_turn/Di;

Massx=Txx*txx/Isp/g0*12
Massy=Tyy*tyy/Isp/g0*12
Massy_turn=Tyy_turn*t_turn/Isp/g0*12


%%%% Delta V orbit

Isp=235;
t_burn=600;
delta_V=19.2;
M0=10000;
g0=9.81;

M_p=M0-M0/(exp(delta_V/(Isp*g0)));
F=Isp*g0*M_p/t_burn
M_p=M_p*2

M_tot=Massx+Massy+M_p

V=M_tot/1.002*1.2
M_tank=2.7086*10^-8 *V^3 -6.1703*10^-5 *V^2 +6.66290*10^-2 *V +1.3192
M_thruster=8*1.6

Mtank_fuel=M_tot+M_tank+M_thruster