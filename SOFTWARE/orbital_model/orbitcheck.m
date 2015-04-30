clc
clear all 
close all

% load constants
constants

do_plot = false;
write_to_file = true;


ry = 10* R_m; %[m]
v = 7000; %[m/s]
dt = 1;

CLCD = 0.3; %[-]
refinement_steps = 40;


% changing variables

rx = -4.14e6:-2e3:-4.2e6;
CD = 1.05:0.05:1.5;


cc = parula(length(CD)+3);
if do_plot
    figure('name','orbits')
end

if write_to_file
    fid = fopen('orbit_true_or_false.txt','a');
end
for k = 1:length(CD)
    CL = CLCD * CD(k);
    % store first value of crashed..
    rx_crashed = rx(1);
    for i=1:length(rx)
        [out] = orbit_selection(rx(i),ry,CD(k),v,dt,CL,R_m,Omega_m,S,m,G,M_mars,h_atm,crash_margin,g_earth);
                            
        if write_to_file
            fprintf(fid,'%11.1f %4.2f %d %d %d %f \n',rx(i),CD(k),out.inatmos,out.crash,out.inorbit,out.maxaccel);
        end
        
        % store last rx value if crashed
        if (out.crash)
            rx_crashed = rx(i);
        end
        
        % check if not crashed and not in orbit
        if (out.crash == false) && (out.inorbit == false)
           
            rx_refine = linspace(rx_crashed,rx(i),refinement_steps+2);
            for j = 2:(length(rx_refine)-1)
                [out] = orbit_selection(rx_refine(j),ry,CD(k),v,dt,CL,R_m,Omega_m,S,m,G,M_mars,h_atm,crash_margin,g_earth);
                if write_to_file
                    fprintf(fid,'%11.1f %4.2f %d %d %d %f \n',rx_refine(j),CD(k),out.inatmos,out.crash,out.inorbit,out.maxaccel);
                end
                
                if do_plot && (out.crash == false)
                    t = 0:dt:(length(out.R)*dt-dt);
                    subplot(3,1,1)
                    hold on
                    grid on
                    plot(t,sqrt(out.R(:,1).^2 + out.R(:,2).^2 + out.R(:,3).^2),'color',cc(k,:));
                    subplot(3,1,2)
                    hold on
                    grid on
                    plot(t,sqrt(out.V(:,1).^2 + out.V(:,2).^2 + out.V(:,3).^2),'color',cc(k,:));
                    subplot(3,1,3)
                    hold on
                    grid on
                    plot(t,sqrt(out.a(:,1).^2 + out.a(:,2).^2 + out.a(:,3).^2),'color',cc(k,:));
                end
            end
            break;
        end
        
        
        
    end
end


if write_to_file
    fclose(fid);
end










% % circle plot:
% theta_plot = 0:0.01:2*pi;
% radius_mars = ones(1,length(theta_plot)) * R_m;
% radius_mars_atmos = ones(1,length(theta_plot)) * (R_m + 150000);
% figure('name','Orbit')
% grid on
% axis equal
% hold on
% plot(R(:,1),R(:,2))
% polar(theta_plot,radius_mars,'r');
% polar(theta_plot,radius_mars_atmos,'g')
% 
% 
% figure('name','parameters over time')
% subplot(3,1,1)
% Rm = sqrt(R(:,1).^2 + R(:,3).^2 + R(:,2).^2);
% plot(t,Rm)
% grid on
% subplot(3,1,2)
% Vm = sqrt(V(:,1).^2 + V(:,3).^2 + V(:,2).^2);
% plot(t,Vm)
% grid on
% subplot(3,1,3)
% am = sqrt(a(:,1).^2 + a(:,3).^2 + a(:,2).^2);
% plot(t,am)
% grid on