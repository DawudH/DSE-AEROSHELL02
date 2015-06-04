function [ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x )
%OPTIMIZATIONWRAPPER Wrapper for the assessGeometry function

LoverD = -0.3;
q = 40;

radius = x(3);
height = x(2);

skewness = 0;%x(1);
poly = [x(4:end),0,0];
[ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( skewness, height, radius, poly, q, LoverD );

% disp('Score:strcat(num2str(score))

end

