function [ failure, mod ] = constraintsWrapper( x )
%constraintsWrapper Wrapper for the constraints function
params = globalParams();

if params.allowskewness
    skewness = x(1);
else
    skewness = 0;
end

heightfactor = x(2);
height = params.radius * heightfactor;

poly = [x(3:end),0,0];

[ failure, mod ] = assessGeometry( skewness, height, params.radius, poly, params.q, params.LoverD );

% disp('Score:strcat(num2str(score))


end