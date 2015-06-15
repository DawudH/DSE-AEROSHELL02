clear; clc;
% load('aeroshapes\iteration1_0.mat');
% clearvars -except x;
saveflag = true;
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
b(3) = -params.minheightfactor;

%heightfactor < maxheightfactor
A(4,2) = 1;
b(4) = params.maxheightfactor;

opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv});
% opts = gaoptimset('PlotFcns',{@gaplotstopping, @gaplotpareto});
opts = gaoptimset(opts, 'UseParallel', true);
opts = gaoptimset(opts, 'PopulationSize', 200);

% opts = gaoptimset(opts, 'InitialPopulation', x);
[x,Fval,exitFlag,Output, population, scores] = ga(@optimizationWrapper,xlength,A,b,[],[],[],[],[],opts);
% [x,Fval,exitFlag,Output, population, scores] = gamultiobj(@optimizationWrapper,xlength,A,b,[],[],[],[],[],opts);

[ score, mod, CoGshift, CD, failed ] = optimizationWrapper( x );
mod.plotCp(true, false);
if saveflag
    save('aeroshapes/optimizationoutputoutput');
end
