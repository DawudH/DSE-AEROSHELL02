clear; clc;

radius = 6;
height = 3;
skewness = 0;
poly = [8.4005    3.8279    4.3011         0         0];
q = 30;
LoverD = -0.3;

[ score, Cmalpha, CDA, penalty, mod  ] = assessGeometry( radius, height, skewness, poly, q, LoverD );
score
penalty
plot(mod.alpha_array, mod.CLCD_array)