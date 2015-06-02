function [ score, Cmalpha, CDA, penalty, mod  ] = assessGeometry( radius, height, skewness, poly, q, LoverD )
%ASSESSGEOMETRY Summary of this function goes here

    
    %Initialise
    alpha0 = 0; %degrees
    dalpha = 1; %degrees
    alphaend = 30; %degrees
    a = 150;
    gamma = 1.29;
    rho = 1e-5;
    T = 150;

    [TriGeom, A] = ParaGeom(q, skewness, radius, height, poly);
    center = [3 0 0];


    [ TriGeom, A, center ] = generategeometry( q );
    geom = aeroGeometry(TriGeom, A);
    mod = modnewtonian(geom, gamma, a, center, rho, T);
    mod = mod.alphasweep(a*40, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

    Cmalpha = 0;
    CDA = 0;
    penalty = 0;
    
    %Calculate trim angle
    helparray = mod.CLCD_array-LoverD;

    if sum(helparray>0)>0 && sum(helparray<0)>0
        [~,alphatrimindex] = min(abs(helparray));
        
        %Calculate performance
        Cmalpha = (mod.CMA_aero_array(2,alphatrimindex+1)-mod.CMA_aero_array(2, alphatrimindex))/dalpha;
        CDA = mod.CRA_aero_array(1,alphatrimindex);
    
    else
        penalty = 1;
    end
        

    
    %Calculate score
    Cmalphafactor = 1;
    CDAfactor = 1;
    penaltyfactor = -1000;
    
    score = (Cmalphafactor * Cmalpha) + (CDAfactor * CDA) + (penaltyfactor * penalty);
end

