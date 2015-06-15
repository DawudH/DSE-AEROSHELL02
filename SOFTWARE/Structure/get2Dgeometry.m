function [xshape, angle] = get2Dgeometry(z)
    % v0.x
    x = [1.6775, 0.6102, 1.6195, -1.0254, 2.3241, -0.0920, 0.1023, 0.6761, -0.1451, 0.7528, -0.4022, 2.0092, 2.7850, 2.9518];

    heightfactor = x(2);
    radius = 6;
    height = radius * heightfactor;

    poly = [x(3:end),0,0];    

    deltaz = 0.00001;
    xshape = polyval(poly, z/radius)/sum(poly)*height;
    
    dx = polyval(poly, (z+deltaz)/radius)/sum(poly)*height;
    
    
    angle = atand((dx-xshape)/deltaz);

    
end
    