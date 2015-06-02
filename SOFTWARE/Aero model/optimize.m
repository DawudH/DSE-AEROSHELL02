clear; clc;
order = 9;
x0 = ones(order,1);
A = diag(x0);
b = inf * x0;

% fmincon(@optimizationWrapper, x0, A, b)


opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv}, 'UseParallel', true);
[x,Fval,exitFlag,Output, population, scores] = ga(@optimizationWrapper,order,[],[],[],[],[],[],[],opts);

[ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x );
save('outputfiles/optimizationoutput');