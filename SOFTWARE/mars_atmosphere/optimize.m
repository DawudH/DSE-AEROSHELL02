% opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping, @gaplotbestindiv});
opts = gaoptimset('PlotFcns',{@gaplotstopping, @gaplotpareto});
opts = gaoptimset(opts, 'UseParallel', true);
opts = gaoptimset(opts, 'PopulationSize', 48);
A = diag([-1 -1 -1]);
b = [-0.5; -500; -3];
% [x,Fval,exitFlag,Output, population, scores] = ga(@getmass,3,A,b,[],[],[],[],[],opts);
[x,Fval,exitFlag,Output, population, scores] = gamultiobj(@getmass,3,A,b,[],[],[],[],[],opts);
