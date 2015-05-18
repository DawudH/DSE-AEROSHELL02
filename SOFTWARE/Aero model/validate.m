% close all;

cases.deg30cone10alpha = 'deg30cone10alpha';
cases.deg30cone20alpha = 'deg30cone20alpha';
cases.sphere = 'sphere12m';
cases.deg60cone = 'deg60cone';
cases.apollo = 'apollo';
cases.irve = 'irve';

shapecase = cases.apollo;


switch shapecase
    
    case cases.irve
        clc; close all; clear;
        qmax_measured = 14.4;
        a = 329.799;
        gamma = 1.4;
        center = zeros(3,1);
        rho = 0.000977525;
        T = 270.650;
        q = 15;
        M = 7;
        [ coords, tri, A ] = generategeometry( 'irvevalidation', q );

        mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        mod = mod.calcAeroangle(M*a,deg2rad(0),0);
        mod.plotCp(true, false);
        disp('Validation of maximum heat flux of IRVE-3');
        disp('Validation data provided Dillman IRVE-3 flight performance');
        disp('Settings: ');
        disp(strcat('gamma: ', num2str(gamma)));
        disp(strcat('M: ', num2str(M)));
        disp(strcat('q_max as Newtonian: ', num2str(mod.qmax_array(end))));
        disp(strcat('q_max measured: ', num2str(qmax_measured)));
        disp(strcat('q_newton/q_measured: ', num2str(mod.qmax_array(end)/qmax_measured)));
        
    case cases.deg60cone
        clc; close all; clear;
        CDmeasured = 0.58;
        a = 300;
        gamma = 1.4;
        center = zeros(3,1);
        rho = 1e-3;
        T = 150;
        q = 100;
        M = 8;
        [ coords, tri, A ] = generategeometry( 'deg60cone', q );

        mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        mod = mod.calcAeroangle(M*a,deg2rad(0),0);
        
        disp('Validation of CD of a 60 degree (30deg halfangle) cone');
        disp('Validation data provided by Anderson fundamentals page 783 and Stevens, 1950');
        disp('Settings: ');
        disp(strcat('gamma: ', num2str(gamma)));
        disp(strcat('M: ', num2str(M)));
        disp(strcat('CD as Newtonian: ', num2str(mod.CR_aero_array(1,end))));
        disp(strcat('CD as measured: ', num2str(CDmeasured)));
        disp(strcat('CDnewton/CDmeasured: ', num2str(mod.CR_aero_array(1,end)/CDmeasured)));
        
    case cases.sphere
        clc; close all; clear;
        CDmeasured = 0.92;
        a = 300;
        gamma = 1.4;
        center = zeros(3,1);
        rho = 1e-3;
        T = 150;
        q = 80;
        M = 10;
        [ coords, tri, A ] = generategeometry( 'sphere12m', q );

        mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        mod = mod.calcAeroangle(M*a,deg2rad(0),0);
        
        disp('Validation of CD of a sphere of radius');
        disp('Validation data provided by Anderson fundamentals page 783 and "Sphere drag measurements ...", 1966');
        disp('Settings: ');
        disp(strcat('gamma: ', num2str(gamma)));
        disp(strcat('M: ', num2str(M)));
        disp(strcat('CD as Newtonian: ', num2str(mod.CR_aero_array(1,end))));
        disp(strcat('CD as measured: ', num2str(CDmeasured)));
        disp(strcat('CDnewton/CDmeasured: ', num2str(mod.CR_aero_array(1,end)/CDmeasured)));

    case cases.apollo
        clear; clc; close all;
        a = 300;
        M = 10.18; % As in the paper
        V = a*14.9;
        gamma = 1.67;
        center = zeros(3,1); % not important
        rho = 1e-3; %not important
        T = 150; % not important
        q = 21;
        CP_apollo = open('CP_apollo.mat');
        zvalidation = CP_apollo.CP_apollo(:,1);
        Cpvalidation = CP_apollo.CP_apollo(:,2);
        
        [ coords, tri, A ] = generategeometry( 'apollovalidation', q );

        modn = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        modn = modn.calcAeroangle(V,deg2rad(0),0);
        modn.plotCp(true, false);
        zsample = 0; % Just one sample
        points = modn.getPointsOnXYPlane(zsample);
        yvalues = modn.coords(2,points);
        Cps = zeros(size(points));
        for i = 1:length(Cps)
            Cps(i) = modn.calcCpOnPoint(points(i));
        end
        
        
        
        modnewtonianmatrix = [yvalues',Cps', points'];
        modnewtonianmatrix = sortrows(modnewtonianmatrix, 1);
%         Rb = 1.956;
%         distancevalues = zeros(size(points'));
%         for i = length(points)/2+1.5:length(points)
%             distancevalues(i) = norm(modn.coords(:,modnewtonianmatrix(i-1,3))-modn.coords(:,modnewtonianmatrix(i,3)))+distancevalues(i-1);
%         end
%         for i = length(points)/2-0.5:-1:1
%             distancevalues(i) = -norm(modn.coords(:,modnewtonianmatrix(i+1,3))-modn.coords(:,modnewtonianmatrix(i,3)))+distancevalues(i+1);
%         end
        
        Cpvalidation = Cpvalidation * modn.Cpmax_array(end);
        

        disp('Validation of Cp for Apollo');
        disp('Validation data provided by Bertin, page 292');
        disp('Settings:');
        disp(strcat('M: ', num2str(M)));
        disp(strcat('gamma: ', num2str(gamma)));
        
        figure;
        plot(modnewtonianmatrix(:,1)*1.095/max(modnewtonianmatrix(:,1)), modnewtonianmatrix(:,2));
        hold on;
        plot(zvalidation, Cpvalidation, '*');
        xlabel('y (m)');
        ylabel('Cp');
        legend('Modified Newtonian', 'Measured');     
        
    case cases.deg30cone10alpha % See the work by Cleary, JW in mendeley, and Bertin page 287
        clear; clc; close all;
        a = 300;
        M = 14.9; % As in the paper
        V = a*14.9;
        gamma = 1.67;
        center = zeros(3,1); % not important
        rho = 1e-3; %not important
        T = 150; % not important
        q = 80;
        
        betavalidation = 0:30:180;
        Cpvalidation = [0.359,0.326,0.234,0.112,0.06,0.037,0.039];
        
        [ coords, tri, A ] = generategeometry( 'deg30cone', q );

        modn = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        modn = modn.calcAeroangle(V,deg2rad(10),0);
        modn.plotCp(true, false);
        xsample = -3.968256554883363; % Just one sample
        points = modn.getPointsOnYZPlane(xsample);
        beta_points = mod(atan2(modn.coords(2,points), modn.coords(3,points)),2*pi);
        beta_points = rad2deg(beta_points(1:end/2+1));
        Cps = zeros(size(beta_points));
        for i = 1:length(Cps)
            Cps(i) = modn.calcCpOnPoint(points(i));
        end

        disp('Validation of Cp for a 30-degree cone');
        disp('Validation data provided by Bertin, page 287');
        disp('Settings:');
        disp(strcat('M: ', num2str(M)));
        disp(strcat('gamma: ', num2str(gamma)));
        
        figure;
        plot(beta_points, Cps);
        hold on;
        plot(betavalidation, Cpvalidation, '*');
        xlabel('beta (deg)');
        ylabel('Cp');
        legend('Modified Newtonian', 'Measured');     
        
    case cases.deg30cone20alpha % See the work by Cleary, JW in mendeley, and Bertin page 287
        clear; clc; close all;
        a = 300;
        M = 14.9; % As in the paper
        V = a*14.9;
        gamma = 1.67;
        center = zeros(3,1); % not important
        rho = 1e-3; %not important
        T = 150; % not important
        q = 80;
        
        betavalidation = 0:30:180;
        Cpvalidation = [0.63, 0.550, 0.327, 0.093, 0.028, 0.009, 0.014];
        
        [ coords, tri, A ] = generategeometry( 'deg30cone', q );

        modn = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
        modn = modn.calcAeroangle(V,deg2rad(20),0);
        modn.plotCp(true, false);
        xsample = -3.968256554883363; % Just one sample
        points = modn.getPointsOnYZPlane(xsample);
        beta_points = mod(atan2(modn.coords(2,points), modn.coords(3,points)),2*pi);
        beta_points = rad2deg(beta_points(1:end/2+1));
        Cps = zeros(size(beta_points));
        for i = 1:length(Cps)
            Cps(i) = modn.calcCpOnPoint(points(i));
        end

        disp('Validation of Cp for a 30-degree cone');
        disp('Validation data provided by Bertin, page 287');
        disp('Settings:');
        disp(strcat('M: ', num2str(M)));
        disp(strcat('gamma: ', num2str(gamma)));
        
        figure;
        plot(beta_points, Cps);
        hold on;
        plot(betavalidation, Cpvalidation, '*');
        xlabel('beta (deg)');
        ylabel('Cp');
        legend('Modified Newtonian', 'Measured');            
    otherwise
end