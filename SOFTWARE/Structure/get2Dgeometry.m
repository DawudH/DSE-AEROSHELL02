function [xshape, angle] = get2Dgeometry(z)
    % v0.x
%     x = [1.6775, 0.6102, 1.6195, -1.0254, 2.3241, -0.0920, 0.1023, 0.6761, -0.1451, 0.7528, -0.4022, 2.0092, 2.7850, 2.9518];
    % v1.0
%     x = [0.7233    0.3478   -0.1451    0.8742    0.8407   -0.5954    0.6158   -2.5585    0.5504    1.1762   -0.9721   -0.3460   0.9269    0.6970];
    % V1.1
    x = [0.9769    0.3166   -0.4273    1.0013   -0.0052    0.1828   -0.8995   -1.5853    2.0548    1.1669];
    
    skewness = x(1);
    heightfactor = x(2);
    radius = 6;
    height = radius * heightfactor;

    poly = [x(3:end),0,0];    

    deltaz = 0.00001;
%     xshape1 = polyval(poly, (z-skewness)/radius)/sum(poly)*height;
%     xshape2 = polyval(poly, (z+skewness)/radius)/sum(poly)*height;
%     xshape = [xshape([
    
    dx = polyval(poly, (z+deltaz)/radius)/sum(poly)*height;
    
    
    angle = atand((dx-xshape)/deltaz);

    
end
    