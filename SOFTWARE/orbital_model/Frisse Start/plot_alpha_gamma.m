clear 
close all
clc

% open the file
load('alpha_gamma_final.mat');


variables


% get curve for the flyby limit
    % get flyby points
    index = find(results.flyby == true);
    temp_gamma = results.GAMMA(index);
    temp_alpha = results.ALPHA(index);
    different_alpha = unique( temp_alpha );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_alpha)
        
        % get location of the alpha
        temp_index = find(temp_alpha == different_alpha(i));
        
        flyby.gamma(i,1) = max(temp_gamma(temp_index));
        flyby.alpha(i,1) = different_alpha(i);
        
        
    end
    
    % fit curve
    flyby.polyfit = polyfit(flyby.gamma,flyby.alpha,2);
    
    
% get curve for the accel limit
    % get accel points
    index = find(results.A > g_earth*3);
    temp_gamma = results.GAMMA(index);
    temp_alpha = results.ALPHA(index);
    different_alpha = unique( temp_alpha );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_alpha)
        
        % get location of the alpha
        temp_index = find(temp_alpha == different_alpha(i));
        
        accel.gamma(i,1) = min(temp_gamma(temp_index));
        accel.alpha(i,1) = different_alpha(i);
        
        
    end
    
    % fit curve
    accel.polyfit = polyfit(accel.gamma,accel.alpha,2);
    
    
% get curve for the crash limit
    % get crash points
    index = find(results.crash == true);
    temp_gamma = results.GAMMA(index);
    temp_alpha = results.ALPHA(index);
    different_alpha = unique( temp_alpha );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_alpha)
        
        % get location of the alpha
        temp_index = find(temp_alpha == different_alpha(i));
        
        crash.gamma(i,1) = min(temp_gamma(temp_index));
        crash.alpha(i,1) = different_alpha(i);
        
        
    end
    
    % fit curve
    crash.polyfit = polyfit(crash.gamma,crash.alpha,2);

    
    
    
%% make plot
figure('name','fitted to datapoints')
hold on
plot(flyby.gamma,polyval(flyby.polyfit,flyby.gamma)*180/pi)
plot(accel.gamma,polyval(accel.polyfit,accel.gamma)*180/pi)
plot(crash.gamma,polyval(crash.polyfit,crash.gamma)*180/pi)
grid on

% extrapolated
cc = parula(4);
figure('name','fitted to datapoints and extrapolated')
hold on
%plot([flyby.gamma; accel.gamma; crash.gamma],[flyby.alpha; accel.alpha; crash.alpha]*180/pi,'o')
plot(21.64:0.025:22,polyval(flyby.polyfit,21.64:0.025:22)*180/pi,'o-','color',cc(1,:))
plot(21.64:0.025:22,polyval(accel.polyfit,21.64:0.025:22)*180/pi,'*-','color',cc(2,:))
plot(21.64:0.025:21.93,polyval(crash.polyfit,21.64:0.025:21.93)*180/pi,'d-','color',cc(3,:))
grid on
ylabel('$\alpha$  $\left[deg\right]$','interpreter','latex')
xlabel('$\gamma$  $\left[deg\right]$','interpreter','latex')
legend('Flyby limit','Acceleration limit','Crash limit')
phi_0_CL_15.crash.polyfit = crash.polyfit;
phi_0_CL_15.accel.polyfit = accel.polyfit;
phi_0_CL_15.flyby.polyfit = flyby.polyfit;
save('phi_0_CL_15.mat','phi_0_CL_15')
matlab2tikz('.\LaTeX\alpha_gamma.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
