variables
param = hk_Parameters(R,V,G,M_mars);
theta = -pi:pi/100:pi;
a = param.a;
e = param.e;
r = a*(1-e^2)./(1+e*cos(theta));
figure()
plot(r,theta)