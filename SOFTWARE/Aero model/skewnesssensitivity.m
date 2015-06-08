load('aeroshapes/sombrero_maxCMalpha_order12_fixedheight.mat');
clearvars -except x mod;

params = globalParams();
 
heightfactor = x(2);
poly = x(3:end);
poly = mod.geom.poly;


skewnessarray = 0:0.5:3;

modarray = {};

for i = 1:length(skewnessarray);
    [ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( skewnessarray(i), heightfactor, params.radius, poly, params.q, params.LoverD );
    modarray{i} = mod;
end

figure;
hold on;
for i = 1:7
    mod = modarray{i};
    plot(mod.CLCD_array(1:end-2), mod.CG_offset(1:end-2));
end

