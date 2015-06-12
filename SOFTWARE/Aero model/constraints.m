function [ failure, mod ] = constraints( skewness, height, radius, poly, q, LoverD )
%ASSESSGEOMETRY Assess a geometry for it's performance

    params = globalParams();

    %% Initialise
    
    % Angle of attack values
    alpha0 = 0; %degrees
    dalpha = 2; %degrees
    alphaend = alphamax; %degrees
    beta = 0;
    phi = 0;
    
    % Mission variables
    a = 150;
    M = 40;
    V = a*M;
    gamma = 1.29;
    rho = 1e-5;
    T = 150;
    mod = NaN;
    

    %% Shape failure criteria
    
    % Initialise objective functions
    failed = false;    
    
    % derivative is bigger than zero everywhere
    xtest = 0:0.01:1;
    if sum(polyval(polyder(poly), xtest)<0)>0
        warning('Derivative criteria failure');
        failed = true;
    end
    
    % capsule fits in the outer body
    if params.radius - skewness <= 2.5
        warning('radius - skewness < 2.5');
        failed = true;
    end
    
    if skewness < 0
        warning('skewness < 0');
        failed = true;
    end
    
    % Height is larger than 0
    if height <= 0
        warning('height <= 0');
        failed = true;
    end
    
    if height > 3*radius
        warning('height > 3*radius');
        failed = true;
    end    
    
    
    %% If not failed
    if ~failed
%         CoGheight = 3;
%         r_capsule = params.r_capsule;
%         skewnessz = skewness * r_capsule/params.radius;
%         xcog = CoGheight + max(polyval(poly, (r_capsule+skewnessz)/params.radius)/sum(poly)*height, polyval(poly, (r_capsule-skewnessz)/radius)/sum(poly)*height);
%         center = [xcog, 0, 0];
        center = [0, 0, 0];

        % Calculate aerodynamic properties
        [TriGeom, A] = ParaGeom(q, skewness, radius, height, poly);
        geom = aeroGeometry(TriGeom, A);
        mod = modnewtonian(geom, gamma, a, center, rho, T);
        mod = mod.alphasweep(V, beta, phi, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

        %% Assess performance failure criteria

        %CLCD is achieved
        helparray = mod.CLCD_array-LoverD;
        if sum(helparray>0)==0 || sum(helparray<0)==0
            warning('L over D criteria failure');
            failed = true;
        end
                
    end
    %% Calculate score
    
    if failed
        failure = 100;
        x = [skewness, height/params.radius, poly(1:end-2)]
    else
        failure = -100;
    end
end