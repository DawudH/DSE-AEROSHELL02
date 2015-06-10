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

% to what limit to search?
% only one of them can be true
flybylimit = false;
orbitlimit = true;
accellimit = false;

alpha = -5:2.5:25;
gamma_range = [21, 22];


gamma_accuracy = 0.0005;
load('alpha_gamma_phi_0_CL_15.mat');
k = length(results.GAMMA) + 1;

for j = 1:length(alpha)
    
    notdone = true;
    n = 1;
    
    
    while notdone
        
        control.alpha_init = alpha(j)*pi/180; % rad
        
        if n <= 2
            gamma = gamma_range(n);
        end
                
        %%function
        [out] = full_orbit(R, V, A, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits, use_alpha_profile,r,v,theta0,gamma,hypkep,Crho,control.alpha_init,control.dalphadt);
        
        results.ALPHA(k) = control.alpha_init;
        results.GAMMA(k) = gamma;
        results.A(k) = norm(max(out.A_aero));
        results.orbit(k) = out.c.orbit;
        results.flyby(k) = out.c.flyby;
        results.crash(k) = out.c.crash;
                
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
save('alpha_gamma_phi_0_CL_15.mat','results')
figure(1)
hold on
grid on
xlim([21.8 22.1])
plot(results.GAMMA(find(results.flyby == true)),results.ALPHA(find(results.flyby == true))*180/pi,'o','color','c')
plot(results.GAMMA(find(results.orbit == true)),results.ALPHA(find(results.orbit == true))*180/pi,'x','color','b')
plot(results.GAMMA(find(results.crash == true)),results.ALPHA(find(results.crash == true))*180/pi,'x','color','g')
plot(results.GAMMA(find(results.A > g_earth*3)),results.ALPHA(find(results.A > g_earth*3))*180/pi,'o','color','r')


    
