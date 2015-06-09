clear 
close all
clc

% open the file
load('alpha_gamma_phi_60.mat');


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
figure('name','fitted to datapoints and extrapolated')
hold on
%plot([flyby.gamma; accel.gamma; crash.gamma],[flyby.alpha; accel.alpha; crash.alpha]*180/pi,'o')
plot(21.74:0.001:22,polyval(flyby.polyfit,21.74:0.001:22)*180/pi)
plot(21.85:0.001:22,polyval(accel.polyfit,21.85:0.001:22)*180/pi)
plot(21.74:0.001:21.93,polyval(crash.polyfit,21.74:0.001:21.93)*180/pi)
grid on

phi_60.crash.polyfit = crash.polyfit;
phi_60.accel.polyfit = accel.polyfit;
phi_60.flyby.polyfit = flyby.polyfit;
save('phi_60.mat','phi_60')

