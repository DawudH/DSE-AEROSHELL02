
z = 0.01:0.05:6;

[xshape,angle] = get2Dgeometry(z);
figure(7)
plot(z,xshape-0.25)
xlim([2.5 6])

theta = atan(xshape./z)*180/pi;