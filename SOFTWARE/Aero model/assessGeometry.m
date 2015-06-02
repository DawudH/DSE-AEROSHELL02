function [ score, Cmalpha, CDA, failed, mod  ] = assessGeometry( radius, height, skewness, poly, q, LoverD )
%ASSESSGEOMETRY Assess a geometry for it's performance

    
    %% Initialise
    
    % Angle of attack values
    alpha0 = 0; %degrees
    dalpha = 5; %degrees
    alphaend = 30; %degrees
    beta = 0;
    
    % Mission variables
    a = 150;
    M = 40;
    V = a*M;
    gamma = 1.29;
    rho = 1e-5;
    T = 150;
    
    % Calculate CoG height
    CoGheight = 3;
    r_capsule = 2.5;
    skewnessz = skewness * r_capsule/radius;
    xcog = CoGheight + max(polyval(poly, (r_capsule+skewnessz)/radius)/sum(poly)*height, polyval(poly, (r_capsule-skewnessz)/radius)/sum(poly)*height);
    center = [xcog, 0, 0];

    % Calculate aerodynamic properties
    [TriGeom, A] = ParaGeom(q, skewness, radius, height, poly);
    geom = aeroGeometry(TriGeom, A);
    mod = modnewtonian(geom, gamma, a, center, rho, T);
    mod = mod.alphasweep(V, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

    
    % Initialise objective functions
    Cmalpha = 0;
    CDA = 0;
    failed = false;
    
    
    %% Assess failure criteria
    % derivative is bigger than zero everywhere
    xtest = 0:0.01:1;
    if sum(polyval(polyder(poly), xtest)<0)>0
        warning('Derivative criteria failure');
        failed = true;
    end
    
    % capsule fits in the outer body
    if radius - skewness < 2.5
        warning('Capsule radius criteria failure');
        failed = true;
    end
    
    %CLCD is achieved
    helparray = mod.CLCD_array-LoverD;
    if sum(helparray>0)==0 || sum(helparray<0)==0
        warning('L over D criteria failure');
        failed = true;
    end
    
    %% Calculate trim angle and function values iff no failure criterion was met
    if ~failed
        
        for i = 1:length(helparray)-1
            if mod.CLCD_array(i+1)<LoverD
                alphatrimindex = i;
                break;
            end
        end
%         [~,alphatrimindex] = min(abs(helparray));
%         alphatrimindex = alphatrimindex(1);
        dCLCDdalpha = (mod.CLCD_array(alphatrimindex+1)-mod.CLCD_array(alphatrimindex))/(mod.alpha_array(alphatrimindex+1)-mod.alpha_array(alphatrimindex));
        realalphatrim = (LoverD-mod.CLCD_array(alphatrimindex))/dCLCDdalpha + mod.alpha_array(alphatrimindex);

        mod.calcAeroangle(V, realalphatrim, beta);
        mod.calcAeroangle(V, realalphatrim+0.001, beta);

        %Calculate performance
        Cmalpha = (mod.CMA_aero_array(2,end)-mod.CMA_aero_array(2, end-1))/(mod.alpha_array(end)-mod.alpha_array(end-1));
        CDA = mod.CRA_aero_array(1,end-1);
    end
        
    %% Calculate score
    Cmalphafactor = 0.01;
    CDAfactor = -1;
    penaltyfactor = +1e10;
    
    score = (Cmalphafactor * Cmalpha) + (CDAfactor * CDA) + (penaltyfactor * failed);
end