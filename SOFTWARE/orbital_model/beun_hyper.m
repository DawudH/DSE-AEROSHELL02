a = 100;
e = 3;
theta = -pi:pi/100:pi;
r = a*(1-e^2)./(1+e*cos(theta));
figure()
plot(r,theta);