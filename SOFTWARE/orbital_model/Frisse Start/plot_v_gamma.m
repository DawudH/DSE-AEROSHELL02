clear 
close all
clc

% open the file
load('v_gamma.mat');

variables


% get curve for the flyby limit
    % get flyby points
    index = find(results.flyby == true);
    temp_gamma = results.GAMMA(index);
    temp_v = results.V(index);
    different_v = unique( temp_v );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_v)
        
        % get location of the alpha
        temp_index = find(temp_v == different_v(i));
        
        flyby.gamma(i,1) = max(temp_gamma(temp_index));
        flyby.v(i,1) = different_v(i);
        
        
    end
    
    % fit curve
    flyby.polyfit = polyfit(flyby.gamma,flyby.v,2);
    
    
% get curve for the accel limit
    % get accel points
    index = find(results.A > g_earth*3);
    temp_gamma = results.GAMMA(index);
    temp_v = results.V(index);
    different_v = unique( temp_v );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_v)
        
        % get location of the alpha
        temp_index = find(temp_v == different_v(i));
        
        accel.gamma(i,1) = min(temp_gamma(temp_index));
        accel.v(i,1) = different_v(i);
        
        
    end
    
    % fit curve
    accel.polyfit = polyfit(accel.gamma,accel.v,2);
    
    
% get curve for the crash limit
    % get crash points
    index = find(results.crash == true);
    temp_gamma = results.GAMMA(index);
    temp_v = results.V(index);
    different_v = unique( temp_v );
    
    % get max gamma for each angle of attack
    for i = 1:length(different_v)
        
        % get location of the alpha
        temp_index = find(temp_v == different_v(i));
        
        crash.gamma(i,1) = min(temp_gamma(temp_index));
        crash.v(i,1) = different_v(i);
        
        
    end
    
    % fit curve
    crash.polyfit = polyfit(crash.gamma,crash.v,2);

    
    
    
%% make plot
figure('name','fitted to datapoints')
hold on
plot(flyby.gamma,polyval(flyby.polyfit,flyby.gamma))
plot(accel.gamma,polyval(accel.polyfit,accel.gamma))
plot(crash.gamma,polyval(crash.polyfit,crash.gamma))
grid on

% extrapolated
figure('name','fitted to datapoints and extrapolated')
hold on
%plot([flyby.gamma; accel.gamma; crash.gamma],[flyby.v; accel.v; crash.v],'o')
gamma = 21.5:0.001:22.5;
plot(gamma,polyval(flyby.polyfit,gamma))
plot(gamma,polyval(accel.polyfit,gamma))
plot(gamma,polyval(crash.polyfit,gamma))
grid on
