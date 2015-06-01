function [out] = aero_conrol(state,control,aero_coef,dt)
    % state should contain:
        % state.a, the accelerataion at this instant
        % state.CL, the CL at this instant
        % state.CD, the CD at this instant
    % control should contain:
        % control.CL_range, the range of CL values to be used
        % control.CLCD, the constant CLCD value
        % control.a, try to maintain a steady value of this acceleration
        % [m/s^2]
        % control.CLa, the dCL/dalpha
        % control.dalpha, the steps at which alpha can be changed each
        % timestep
    
        % determine the error
        
        e = state.h - control.h;
        if (control.error == 0)
            control.error = e;
        end
        error_I = control.error_I + dt * (e + control.error)/2; % trapozoidal integration
        error_D = (e - control.error) / dt; 
        dalpha = ( control.Kp * e * control.dalpha + control.Ki * error_I * control.dalpha + control.Kd * error_D * control.dalpha);
        dalpha_beun = dalpha;

        if (control.dalpha > 0) && (dalpha > control.dalpha)
            dalpha = control.dalpha;
        elseif (control.dalpha > 0) && (dalpha < -control.dalpha)
            dalpha = -control.dalpha;
        elseif (control.dalpha < 0) && (dalpha > -control.dalpha)
             dalpha = -control.dalpha;
        elseif (control.dalpha < 0) && (dalpha < control.dalpha)
             dalpha = control.dalpha;
        end

        
        alpha = state.alpha + dalpha;
        if ( alpha < min(control.alpha_range) )
             alpha = min(control.alpha_range);
        elseif ( alpha > max(control.alpha_range) )
             alpha = max(control.alpha_range);
        end
        
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
        
    
    
    % generate output
    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    out.error = e;
    out.error_I = error_I;
    disp(['error: ' num2str(e) ' error_I: ' num2str(error_I) ' error_D: ' num2str(error_D) ' dalpha_beun: ' num2str(dalpha_beun) ' dalpha: ' num2str(dalpha) ' alpha: ' num2str(alpha*180/pi)])
end