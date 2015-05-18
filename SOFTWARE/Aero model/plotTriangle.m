function [  ] = plotTriangle( mod, n )
%PLOTTRIANGLE Plot the triangle

t = mod.tri(n,:);

plot3(mod.coords(1,t(1)), mod.coords(2,t(1)), mod.coords(3,t(1)), 'o');
plot3(mod.coords(1,t(2)), mod.coords(2,t(2)), mod.coords(3,t(2)), 'o');
plot3(mod.coords(1,t(3)), mod.coords(2,t(3)), mod.coords(3,t(3)), 'o');

end

