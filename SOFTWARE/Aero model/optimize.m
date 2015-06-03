clear; clc;
order = 8;
xlength = order+3;
x0 = ones(xlength,1);

b = zeros(4,1);

% fmincon(@optimizationWrapper, x0, A, b)
A = zeros(4, xlength);
% -radius+skewness<2.5
A(1,1) = 1;
A(1,3) = -1;
b(1) = 2.5;

% skewness>0
A(2,1) = -1;
b(2) = 0;

% height > 0
A(3,2) = -1;
b(3) = 0;

% radius > 0
A(4,3) = -1;
b(4) = 0;

opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv}, 'UseParallel', true);
[x,Fval,exitFlag,Output, population, scores] = ga(@optimizationWrapper,xlength,[],[],[],[],[],[],[],opts);

[ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x );
save('outputfiles/optimizationoutput');