function [ score, mod, CoGshift, CDA, failed ] = optimizationWrapper( x )
%OPTIMIZATIONWRAPPER Wrapper for the assessGeometry function
params = globalParams();

if params.allowskewness
    skewness = x(1);
else
    skewness = 0;
end

heightfactor = x(2);
height = params.radius * heightfactor;

poly = [x(3:end),0,0];

[ score, CoGshift, CDA, failed, mod  ] = assessGeometry( skewness, height, params.radius, poly, params.q, params.LoverD );
score = score(6);
disp(score);
% disp('Score:strcat(num2str(score))

end