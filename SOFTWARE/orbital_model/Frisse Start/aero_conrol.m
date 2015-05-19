function [out] = aero_conrol(state,control,aero_coef)
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
        e = abs(control.a - state.a);
        dalpha = control.Kp * e * control.dalpha;
        if dalpha > control.dalpha
            dalpha = control.dalpha;
        end
        
    if state.a < control.a
    % If state.a < control.a, decrease CL (so decrease alpha) to allow higer acceleration
        
        alpha = state.alpha - dalpha;
        if ( alpha < min(control.alpha_range) )
             alpha = min(control.alpha_range);
        end
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
        
         
    elseif state.a == control.a
        % If state.a = control.a, keep the same CL
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(state.alpha);
    
    else
    % If state.a > control.a, increase CL (so increase alpha) to allow lower acceleration
    
        alpha = state.alpha + dalpha;
        if ( alpha > max(control.alpha_range) )
             alpha = max(control.alpha_range);
        end
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
        
    end
    
    
    % generate output
    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    
end