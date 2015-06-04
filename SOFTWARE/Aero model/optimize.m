clear; clc;

order = 8;
xlength = order+2;

x0 = ones(xlength,1);

b = zeros(2,1);

% fmincon(@optimizationWrapper, x0, A, b)
A = zeros(2, xlength);

% % -radius+skewness<-2.5
% A(1,1) = 1;
% A(1,3) = -1;
% b(1) = -2.5;

% skewness>0
A(1,1) = -1;
b(1) = 0;

% height > 0
A(2,2) = -1;
b(2) = 0;


opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv});
opts = gaoptimset(opts, 'UseParallel', true);
opts = gaoptimset(opts, 'PopulationSize', 50);
[x,Fval,exitFlag,Output, population, scores] = ga(@optimizationWrapper,xlength,A,b,[],[],[],[],[],opts);

[ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x );
save('aeroshapes/optimizationoutput_maxmoment');