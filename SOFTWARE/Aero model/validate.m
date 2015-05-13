% close all;

shapetexts.horizontalplate = 'horizontalplate';
shapetexts.verticalplate = 'verticalplate';
shapetexts.sphere12m = 'sphere12m';
shapetexts.pastille12m15m = 'pastille12m1.5m';
shapetexts.deg60cone = 'deg60cone';
shapetexts.deg30cone = 'deg30cone';

shape = shapetexts.deg30cone;


switch shape
    case shapetexts.deg30cone %See the work by Cleary, JW in mendeley, and Bertin page 287
        a = 300;
        M = 14.9; % As in the paper
        V = a*14.9;
        gamma = 1.67;
        center = zeros(3,1); % not important
        rho = 1e-3; %not important
        T = 150; % not important
        q = 21;
        
        [ coords, tri, A ] = generategeometry( shapetexts.deg30cone, q );

        mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        mod = mod.calcAeroangle(V,deg2rad(10),0);
        mod.plotCp(true, false);
        figure;
        xsample = -7.464101615137754; % Just one sample
        points = mod.getPointsOnYZPlane(xsample);
        beta_points = atan2(mod.coords(3,points), mod.coords(1,points));
        plot(beta_points);
        
%         mod.CR_aero_array
    otherwise
end