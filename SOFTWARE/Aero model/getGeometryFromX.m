function geom = getGeometryFromX(x, q)
    if nargin == 1
        q = 31;
    end
    
    params = globalParams();

    skewness = x(1);


    heightfactor = x(2);
    height = params.radius * heightfactor;

    poly = [x(3:end),0,0];
    
    [TriGeom, A] = ParaGeom(q, skewness, params.radius, height, poly);
    geom = aeroGeometry(TriGeom, A, poly);
    geom.plotGeometry(true, false);
    view(90,0);
end