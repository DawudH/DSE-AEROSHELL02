variables

%%Twan
aVV
aVgamma
aVR
aValpha
agammaV
enz...

%%Sebastiaan
agammagamma = -((V0/R0)-(g0/V0))*sin(gamma0);
agammaR = ((2*g0/V0)-(V0/R0))*cos(gamma0)/R0;
agammaalpha = (cos(sigma0)/(m*V0))*CLalpha*q0*Sref;
agammabeta=-(sin(sigma0)/(m*V0))*CSbeta*q0*Sref;
agammasigma = -(L0/(m*V0))*sin(sigma0);
agammap = 0;
agammaq = 0;
agammar = 0;
aRV = sin(gamma0);
aRgamma = V0*cos(gamma0);
aRR = 0;
aRp = 0;
aRq = 0;
aRr = 0;
aRalpha = 0;
aRbeta = 0;
aRsigma = 0;
apbeta = (1/Ixx)*q0*Sref*bref*CLbeta;
apV = 0;
apgamma = 0;
apR = 0;
app = 0;
apq = 0;
apr = 0;
apalpha = 0;
apsigma = 0;
aqV = (M0/(Iyy*V0))*CmM*q0*Sref*cref;
aqalpha = (1/Iyy)*Cmalpha*q0*Sref*cref;

aalphap = 0;
aalphar = 0;
aalphabeta = 0;
abetaV = (g0/(V0*V0))*cos(gamma0)*sin(sigma0);
abetagamma = (g0/V0)*sin(gamma0)*cos(gamma0);
abetaR = 2*(g0/(R0*V0))*cos(gamma0)*sin(sigma0);
abetap = sin(alpha0);
abetar = -cos(alpha0);
abetabeta = -(1/(m*V0))*CSbeta*q0*Sref;
abetasigma = -(g0/V0)*cos(gamma0)*cos(sigma0);
abetaq = 0;
abetaalpha = 0;

asigmaV = (tan(gamma0)*sin(sigma0)/(m*V0*V0))*(M0*CLM+CL)*q0*Sref;
asigmagamma = (L0/(m*V0))*sin(sigma0);
asigmap = -cos(alpha0);
asigmar = -sin(alpha0);
asigmaalpha = (tan(gamma0)*sin(sigma0)/(m*V0))*CLalpha*q0*Sref;
asigmabeta = (tan(gamma0)*cos(sigma0)/(m*V0))*CSbeta*q0*Sref-L0/(m*V0)+(g0/V0)*cos(gamma0)*cos(sigma0);
asigmasigma = tan(gamma0)*cos(sigma0)*L0/(m*V0);
asigmaR = 0;
asigmaq = 0;