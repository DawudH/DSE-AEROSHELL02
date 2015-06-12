addpath('..\aerodynamic_coefficients');
aerocoef=aeroProperties('irve');

Di=4.5;
Do=12;
L=5;

Ds_cb=3.5;
m_cb=9000;

m_h=1000;
m_inf=200;

Isp=234;
g0=9.81;

mu_dd=5;
mu=45;
alpha_dd=5;
alpha=5;

q=1000;

CM_alphaAL=1.5
Myy_turn=CM_alphaAL*alpha*q;
t_turn=5;

m_cbh=m_h-m_inf;
t_cb=Di-Ds_cb;

Ixx=1/4 * (0.5*m_cb*((Ds_cb+t_cb)^2+(Ds_cb+t_cb)^2) + 0.5*m_cbh*Di^2  +  0.5*m_inf+Do^2);
Iyy=1/12*(m_cb+m_h)*(3/4*((Ds_cb+t_cb)^2+(Ds_cb+t_cb)^2) +L^2);


mu_dd=mu_dd/180*pi;
alpha_dd=alpha_dd/180*pi;
mu=mu/180*pi;
alpha=alpha/180*pi;

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
delta_V=20;
M0=10000;
g0=9.81;

M_p=M0-M0/(exp(delta_V/(Isp*g0)))
F=Isp*g0*M_p/t_burn
M_p=M_p*2

M_tot=Massx+Massy+Massy_turn+M_p
