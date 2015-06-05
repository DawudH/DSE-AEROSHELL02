clear; clc;

params = globalParams();

order = 8;
xlength = order+2;


A = zeros(4, xlength);
b = zeros(4,1);

% skewness<-2.5+r
A(1,1) = 1;
b(1) = -2.5 + params.radius;

% skewness>0
A(2,1) = -1;
b(2) = 0;

% height > 0
A(3,2) = -1;
b(3) = 0;

%heightfactor < 3
A(4,2) = 1;
b(4) = 3;

% opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv, @gaplotpareto});
opts = gaoptimset('PlotFcns',{@gaplotstopping, @gaplotpareto});
opts = gaoptimset(opts, 'UseParallel', true);
opts = gaoptimset(opts, 'PopulationSize', 48);
[x,Fval,exitFlag,Output, population, scores] = gamultiobj(@optimizationWrapper,xlength,A,b,[],[],[],[],[],opts);

[ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x );
mod.geom.plotGeometry(true, false);
save('aeroshapes/optimizationoutput_maxmoment');