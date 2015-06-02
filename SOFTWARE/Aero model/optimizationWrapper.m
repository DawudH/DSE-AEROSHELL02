function [ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x )
%OPTIMIZATIONWRAPPER Wrapper for the assessGeometry function

LoverD = -0.3;
q = 40;

radius = 6; %x(1);
height = 3; %x(2);

skewness = 0; %x(1);
poly = [x(2:end),0,0];
[ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( radius, height, skewness, poly, q, LoverD );

% disp('Score:strcat(num2str(score))

end

