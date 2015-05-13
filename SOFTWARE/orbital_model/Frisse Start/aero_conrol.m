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

    if state.a < control.a
    % If state.a < control.a, decrease CL (so decrease alpha) to allow higer acceleration
    
        alpha = state.alpha - control.dalpha;
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
        
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          CL = state.CL - control.CLa * control.dalpha;
%          % check if CL is out of bounds
%          if ( CL < min(control.CL_range) )
%              CL = min(control.CL_range);
%          end
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    elseif state.a == control.a
        % If state.a = control.a, keep the same CL
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(state.alpha);
    
    else
    % If state.a > control.a, increase CL (so increase alpha) to allow lower acceleration
    
        alpha = state.alpha + control.dalpha;
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
        
%         CL = state.CL + control.CLa * control.dalpha;
%          % check if CL is out of bounds
%          if ( CL > max(control.CL_range) )
%              CL = max(control.CL_range);
%          end
    end
    
    
    % generate output
    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    
end