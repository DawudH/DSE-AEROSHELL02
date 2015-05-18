function [  ] = plotPoint(mod, n )
%PLOTPOINT Plot a point on the 3d plot
plot3(mod.coords(1,n), mod.coords(2,n), mod.coords(3,n), 'o');

end

