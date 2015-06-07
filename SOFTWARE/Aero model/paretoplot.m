for i = 1:length(x)
figure
[ score, mod, Cmalpha, CDA, failed ] = optimizationWrapper( x(i,:) );
mod.geom.plotGeometry(true, false);

end