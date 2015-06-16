clc
clear
close all

%%Input constants & variables
variables

% booleans
use_control = false;
multiple_orbits = false;
use_alpha_profile = false;
export_figures = false;
hypkep = false;

% %% 
% alpha = -10:2.5:20; 
% gamma_range = 21.8:0.01:24;
% k = 1;
% for j = 1:length(alpha)
%     for i = 1:length(gamma_range)
%         
%         control.alpha_init = alpha(j)*pi/180; % rad
%         gamma = gamma_range(i);
%         
%         %%function
%         [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
% 
%         if out.c.orbit
%             results.ALPHA(k) = control.alpha_init;
%             results.GAMMA(k) = gamma;
%             results.A(k) = norm(max(out.A_aero));
%             k = k+1;
%         elseif  out.c.crash
%             break;
%         end
%     end
% end
% 
% save('results.mat','results')
% plot(results.GAMMA,results.ALPHA*180/pi,'x')

%%
flybylimit = false;
orbitlimit = false;
accellimit = true;

gamma_accuracy = 0.005;

load('v_gamma_final.mat');

V_esc = sqrt(G*M_mars * 2 / r);

v_range = 6500:250:7500;
k = length(results.GAMMA) + 1;
k = 1;
for j = 1:length(v_range)
    notdone = true;
    gamma_range = [21.5, 22.5];
    n = 1;
    while notdone
        
        v = v_range(j);
        
        if n <= 2
            gamma = gamma_range(n);
        end
             
        
        
        %%function
        [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);

        results.V(k) = v;
        results.GAMMA(k) = gamma;
        results.A(k) = norm(max(out.A_aero));
        results.orbit(k) = out.c.orbit;
        results.flyby(k) = out.c.flyby;
        results.crash(k) = out.c.crash;

        
        % golden section search (to flyby limit)
        if n >= 2
          
            interval = abs(results.GAMMA(k) - results.GAMMA(k-1)) / 2;
            if interval <= gamma_accuracy
                break;
            end
            
            if (flybylimit)
                if results.flyby(k)  % go to the right..
                    gamma = gamma + interval;
                else % go to the left
                    gamma = gamma - interval;
                end
            end
            
            if (orbitlimit)
                if results.orbit(k) || results.flyby(k)  % go to the right..
                    gamma = gamma + interval;
                else % go to the left
                    gamma = gamma - interval;
                end
                    
            end
            
            if (accellimit)
                if results.A(k) < (3*g_earth) % go to the right..
                    gamma = gamma + interval;
                else % go to the left
                    gamma = gamma - interval;
                end
                    
            end
            disp(num2str(gamma))
                
        end
        n = n + 1;
        k = k+1;
        
    end
end
%%
save('v_gamma_final.mat','results')
figure(1)
hold on
grid on
plot(results.GAMMA(find(results.flyby == true)),results.V(find(results.flyby == true)),'o','color','c')
plot(results.GAMMA(find(results.orbit == true)),results.V(find(results.orbit == true)),'x','color','b')
plot(results.GAMMA(find(results.crash == true)),results.V(find(results.crash == true)),'x','color','g')
plot(results.GAMMA(find(results.A > g_earth*3)),results.V(find(results.A > g_earth*3)),'o','color','r')