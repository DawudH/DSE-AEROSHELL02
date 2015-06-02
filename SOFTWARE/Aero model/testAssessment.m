radius = 6;
height = 3;
skewness = 2;
poly = [1 0 0];
q = 30;
LoverD = -0.3;

[ score, Cmalpha, CDA, penalty, mod  ] = assessGeometry( radius, height, skewness, poly, q, LoverD );
score
plot(mod.alpha_array, mod.CLCD_array)