function error = discretization_check(dt)
% input a range of dt every following dt should be dubble of the previous)

% load constants
variables

% run for different timesteps
R_o = cell(length(dt),1);
error = cell(length(dt),1);
max_error = zeros(length(dt),1);

% plot the orbits
% circle plot:
theta_plot = 0:0.01:2*pi;
radius_mars = ones(1,length(theta_plot)) * R_m;
radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + h_atm);
figure('name','Orbit')
grid on
axis equal
hold on
axis([-(R_m + h_atm)*1.5 (R_m + h_atm)*1.5 -(R_m + h_atm)*1.5 (R_m + h_atm)*1.5])
polar(theta_plot,radius_mars,'r');
polar(theta_plot,radius_mars_atmos,'g')
cc = parula(length(dt)+3);

% Loop through different dt
    for i = 1:length(dt)
        [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt(i), m, Omega_m, S, control, tend, crash_margin, g_earth);
        R_o{i} = out.R;
        len_old = length(R_o{1});
        len_new = length(out.R);
        
        if i > 1
            if (len_old/(2^(i-1)) <= len_new)
                k = 1:2^(i-1):len_old;
                error{i} = out.R(1:length(k),:) - R_o{1}(k,:);
            else
                k = 1:2^(i-1):len_new*2^(i-1);
                error{i} = out.R(1:length(k),:) - R_o{1}(k,:);
            end
            
            % plot orbit
            plot(out.R(1:length(k),1),out.R(1:length(k),2),'color',cc(i,:))
        else
            error{i} = out.R - out.R;
            % plot orbit
            plot(out.R(:,1),out.R(:,2),'color',cc(i,:))
        end
        max_error(i,1) = max(sqrt(error{i}(:,1).^2 + error{i}(:,2).^2 + error{i}(:,3).^2));
        
        
        
    end

    figure('name', 'Discretization check')
    loglog(dt,max_error,'o')
    grid on

end
