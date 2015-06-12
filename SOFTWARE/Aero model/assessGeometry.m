% function [ score, mod ] = assessGeometry( skewness, height, radius, poly, q, LoverD )
function [ score, CoGshift, CD, failed, mod  ] = assessGeometry( skewness, height, radius, poly, q, LoverD )
%ASSESSGEOMETRY Assess a geometry for it's performance
    params = globalParams();
    x = [skewness, height/radius, poly(1:end-2)];

    %% Initialise
    
    % Angle of attack values
    alpha0 = 0; %degrees
    dalpha = 2; %degrees
    alphaend = params.alphamax; %degrees
    beta = 0;
    phi = 0;
    
    % Mission variables
    a = 150;
    M = 40;
    V = a*M;
    mod = NaN;
    

    %% Shape failure criteria
    
    % Initialise objective functions
    Cmalpha = 0;
    CD = 0;
    CLA = 0;
    CmAtrim = 0;
    absoluteLoverD = 0;
    CoGshift = 0;
    failed = false;    
    
    % derivative is bigger than zero everywhere
    xtest = 0:0.01:1;
    if sum(polyval(polyder(poly), xtest)<0)>0
        disp('Warning: Derivative criteria failure');
        failed = true;
    end
    
    % capsule fits in the outer body
    if radius - skewness <= 2.5
        disp('Warning: radius - skewness < 2.5');
        failed = true;
    end
    
    if skewness < 0
        disp('Warning: skewness < 0');
        failed = true;
    end
    
    % Height is larger than 0, smaller than 2*radius
    if height <= -params.minheightfactor*radius
        disp('Warning: height < 0.0001*radius');
        failed = true;
    end
    
    % Height factor greater than max height factor
    if height > params.maxheightfactor*radius
        warning(strcat('height>',num2str(params.maxheightfactor),'*radius'));
        failed = true;
    end    
    
    % Half-cone angle smaller than allowed
    z = 2.5:0.1:radius;
    xshape = polyval(poly, (z+abs(skewness))/radius)/sum(poly)*height;
    dz = z(2:end)-z(1:end-1);
    dx = xshape(2:end)-xshape(1:end-1);
    ourtheta = atand(dx./dz);
    if any(ourtheta < 90-params.thetamin)
        disp('Warning: derivative < thetamin');
        failed = true;
    end
    
    
    
    %% If not failed
    if ~failed


        % Calculate aerodynamic properties
        [mod, center] = generateGeometry(poly, q, skewness, radius, height);
        mod = mod.alphasweep(V, beta, phi, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

        %% Assess performance failure criteria

        %CLCD is achieved
        helparray = abs(mod.CLCD_array)-abs(LoverD);
        if sum(helparray>0)==0 || sum(helparray<0)==0
            disp('Warning: L over D criteria failure in assessgeometry!');
            failed = true;
        end
        
        %Stability is achieved
        if sum(mod.Cmyalpha - params.Cmalpharequired<=0)==0
            disp('Warning: Cmalpha criteria failure in assessgeometry!');
            failed = true;
        end        


        %% Calculate trim angle and function values iff no failure criterion was met
        if ~failed
            alphatrimindex = -1;
            for i = 1:length(helparray)-1
                if (((abs(mod.CLCD_array(i+1))>abs(LoverD)) && (abs(mod.CLCD_array(i))<abs(LoverD))) || ...
                    ((abs(mod.CLCD_array(i))>abs(LoverD)) && (abs(mod.CLCD_array(i+1))<abs(LoverD))))
                    alphatrimindex = i;
                    break;
                end
            end
            if alphatrimindex == -1
                disp('alphatrimindex==-1');
                disp(x);
            end    
%             alphatrimindex = 1;

            dCLCDdalpha = (abs(mod.CLCD_array(alphatrimindex+1)-mod.CLCD_array(alphatrimindex)))/(mod.alpha_array(alphatrimindex+1)-mod.alpha_array(alphatrimindex));
            realalphatrim = (abs(LoverD)-abs(mod.CLCD_array(alphatrimindex)))/dCLCDdalpha + mod.alpha_array(alphatrimindex);
            
%             realalphatrim = 0;
            
            mod = mod.calcAeroangle(V, realalphatrim, beta, phi);
            mod = mod.calcAeroangle(V, realalphatrim+0.001, beta, phi);
            if mod.Cmyalpha(end) < params.Cmalpharequired
                Cmalpha = mod.Cmyalpha(end);
            else
                failed = true;
            end
            %Calculate performance
            CD = mod.CR_aero_array(1,end-1);
            CmAtrim = mod.CMA_aero_array(2,end-1);
            absoluteLoverD = max(abs(mod.CLCD_array));
            absoluteCLA = max(abs(mod.CRA_aero_array(3,:)));
            CoGshift = abs(mod.CG_offset(end-1));
        end
        
    end
    %% Calculate score
    if failed
        disp('failed');
        disp(x)
        score = [1000, -1000, 1000, -1000, -10000, 1000];
    else
        score = [Cmalpha;CD;CmAtrim;absoluteLoverD;absoluteCLA;CoGshift];
    end

end