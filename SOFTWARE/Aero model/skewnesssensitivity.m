load('aeroshapes/vingerhoed_maxLoverD_order12_varheight.mat');
clearvars -except x mod;

params = globalParams();
 
heightfactor = 6;%x(2);
poly = x(3:end);
poly = mod.geom.poly;


% skewnessarray = 0:0.5:3;
skewnessarray = 1.5;
modarray = {};

for i = 1:length(skewnessarray);
    [ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( skewnessarray(i), heightfactor, params.radius, poly, params.q, params.LoverD );
    modarray{i} = mod;
end

% figure;
% hold on;
mod.geom.plotGeometry(true, false)

% for i = 1:7
%     mod = modarray{i};
%     plot(mod.CLCD_array(1:end-2), mod.CG_offset(1:end-2));
% end

