function [ mod, center ] = generateGeometry(poly, q, skewness, radius, height, a, gamma, rho, T)
    %Specify either 0, 1, 5, or 9 variables
switch nargin
    %Specified only polynomial
    case 1
        q = 30;
        skewness = 0;
        radius = 6;
        height = 3;
        a = 150;
        T = 150;
        gamma = 1.29;
        rho = 1e-5;        

    %Specified only geometry
    case 5
        a = 150;
        T = 150;
        gamma = 1.29;
        rho = 1e-5;
        
    case 9

    otherwise
        disp('Not all arguments given')
        q = 30;
        skewness = 0;
        poly = [1 0 0];
        radius = 6;
        height = 3;
end

    params = globalParams();

    [TriGeom, A] = ParaGeom(q, skewness, radius, height, poly);
    geom = aeroGeometry(TriGeom, A, poly);
    
    CoGcapsuleheight = 3;
    skewnessz = skewness * params.r_capsule/radius;
    radiusheight1 = polyval(poly, (params.r_capsule+skewnessz)/radius)/sum(poly)*height;
    radiusheight2 = polyval(poly, (params.r_capsule-skewnessz)/radius)/sum(poly)*height;
    capsuleCoGHeight = CoGcapsuleheight + max(radiusheight1, radiusheight2);        
    xcog = (1000*geom.centroid(1)-9000*capsuleCoGHeight)/10000;
    center = [xcog, 0, 0];
    mod = modnewtonian(geom, gamma, a, center, rho, T);
end