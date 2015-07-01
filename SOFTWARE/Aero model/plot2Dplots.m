load('C:\Users\Joost\Documents\DSE-AEROSHELL02\SOFTWARE\Aero model\aeroshapes\iteration1_1_orbit.mat');
x(1) = 1.5;
x(2) = 0.5;
geom = getGeometryFromX(x, 61);
geom.plotView();