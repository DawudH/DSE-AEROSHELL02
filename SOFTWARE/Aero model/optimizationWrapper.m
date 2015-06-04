function [ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x )
%OPTIMIZATIONWRAPPER Wrapper for the assessGeometry function
radiusglobal = 6;
LoverDglobal = -0.3;
qglobal = 30;

heightfactor = x(2);
height = radiusglobal * heightfactor;

skewness = 0; %x(1);
poly = [x(3:end),0,0];
[ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( skewness, height, radiusglobal, poly, qglobal, LoverDglobal );

% disp('Score:strcat(num2str(score))

end