clc
clear all 
close all
tic
% load constants
constants

do_plot = false;
write_to_file = true;
file_name = 'orbit_selection_range.txt';


ry = 10* R_m; %[m]
v = 7000; %[m/s]
dt = 1;
refinement_steps = 100;


% changing variables
% 
rx = -4161499:-5e2:-4.8e6;
CD = 1.25;
%CLCD = [-0.25:0.1:0 0:0.05:0.3]; %[-]
CLCD = 0.25;
CL = -0.4;

cc = parula(length(CD)+3);
if do_plot
    figure('name','orbits');
end

rx_first = rx(1);
% file string
filestr = cell(length(v),1);
parfor u = 1:length(v)
    for l = 1:length(CL)
    for k = 1:length(CD)

        % store first value of crashed..
        rx_crashed = rx_first;
        for i=1:length(rx)
            [out] = orbit_selection(rx(i),ry,CD(k),v(u),dt,CL(l),R_m,Omega_m,S,m,G,M_mars,h_atm,crash_margin,g_earth);

            if write_to_file
                filestr{u} = [filestr{u}  sprintf(' %11.1f %4.2f %4.2f %6.1f %d %d %d %f \n',rx(i),CD(k),CL(l),v(u),out.inatmos,out.crash,out.inorbit,out.maxaccel)];
            end

            % store last rx value if crashed
            if (out.crash) || (out.inorbit)
                rx_crashed = rx(i);
            end

            % check if not crashed and not in orbit
            if (out.crash == false) && (out.inorbit == false)
                previous_inorbit = false;
                rx_refine = linspace(rx_crashed,rx(i),refinement_steps+2);
                for j = 2:(length(rx_refine)-1)
                    [out] = orbit_selection(rx_refine(j),ry,CD(k),v(u),dt,CL(l),R_m,Omega_m,S,m,G,M_mars,h_atm,crash_margin,g_earth);
                    if write_to_file
                        filestr{u} = [filestr{u} sprintf(' %f %f %f %f %d %d %d %f \n',rx(j),CD(k),CL(l),v(u),out.inatmos,out.crash,out.inorbit,out.maxaccel)];
                    end

                    if previous_inorbit && (out.inorbit == 0)
                        break;
                    end

                    previous_inorbit = out.inorbit;

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
    end
end
if write_to_file
    fid = fopen(file_name,'a');
    for k = 1:length(v)
        fprintf(fid,'%s',filestr{k});
    end
    fclose(fid);
end

toc






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