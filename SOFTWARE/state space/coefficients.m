% Load variables
[go,m,bref,Sref,V0,M0,rho0,alpha0,gamma0,gammadot0,q0,CL0A,CD0A,CM0A,D0,L0,MY0,R0,sigma0,Ixx,Iyy,Izz]=variables()

aVV = -(2*D0)/(m*V0);
aVgamma = -g0*cos(gamma0);
aVR = 2*g0/R0*sin(gamma0);
aValpha = -1/m*CDalpha*q0*Sref;
agammaV = (-gammadot0+2*V0/R0*cos(gamma0))/V0 + cos(sigma0)*2*L0/(m*V0^2);
aVp = 0;
aVq = 0;
aVr = 0;
aVbeta = 0;
aVsigma = 0;
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
aqgamma = 0;
aqR = 0;
aqp = 0;
aqq = 0;
aqr = 0;
aqbeta = 0;
aqsigma = 0;
arbeta = CNbeta*q0*Sref*bref;
arV = 0;
argamma = 0;
arR = 0;
arp = 0;
arq = 0;
arr = 0;
aralpha = 0;
arsigma = 0;
aalphaV = -g0/V0^2*cos(gamma0)*cos(sigma0)-L0/(m*V0^2);
aalphagamma = -g0/V0*sin(gamma0)*cos(sigma0);
aalphaR = -2*g0/(R0*V0)*cos(gamma0)*cos(sigma0);
aalphaq = 1;
aalphaalpha = 1/(m*V0)*CLalpha*q0*Sref;
aalphasigma = -g0/V0*cos(gamma0)*sin(sigma0);
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


A = [aVV, aVgamma, aVR, aVp, aVq, aVr, aValpha, aVbeta, aVsigma;...
     agammaV, agammagamma, agammaR, agammap, agammaq, agammar, agammaalpha, agammabeta, agammasigma;...
     aRV, aRgamma, aRR, aRp, aRq, aRr, aRalpha, aRbeta, aRsigma;...
     apV, apgamma, apR, app, apq, apr, apalpha, apbeta, apsigma;...
     aqV, aqgamma, aqR, aqp, aqq, aqr, aqalpha, aqbeta, aqsigma;...
     arV, argamma, arR, arp, arq, arr, aralpha, arbeta, arsigma;...
     aalphaV, aalphagamma, aalphaR, aalphap, aalphaq, aalphar, aalphaalpha, aalphabeta, aalphasigma;...
     abetaV, abetagamma, abetaR, abetap, abetaq, abetar, abetaalpha, abetabeta, abetasigma;...
     asigmaV, asigmagamma, asigmaR, asigmap, asigmaq, asigmar, asigmaalpha, asigmabeta, asigmasigma];