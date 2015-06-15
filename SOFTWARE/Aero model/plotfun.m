function state = plotfun(options, state, flag)
    x = state.Population(1,:);
    
    
    params = globalParams();

    skewness = x(1);


    heightfactor = x(2);
    height = params.radius * heightfactor;

    poly = [x(3:end),0,0];
    
    [TriGeom, A] = ParaGeom(31, skewness, params.radius, height, poly);
    geom = aeroGeometry(TriGeom, A, poly);
    geom.plotGeometry(true, false);
end