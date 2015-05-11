function [out] = aero_conrol(state,control)
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
    % If state.a < control.a, decrease CL to allow higer acceleration
         CL = state.CL - control.CLa * control.dalpha;
         % check if CL is out of bounds
         if ( CL < min(control.CL_range) )
             CL = min(control.CL_range);
         end
         
    elseif state.a == control.a
    % If state.a = control.a, keep the same CL
        CL = state.CL;
    
    else
    % If state.a > control.a, increase CL to allow lower acceleration
        CL = state.CL + control.CLa * control.dalpha;
         % check if CL is out of bounds
         if ( CL > max(control.CL_range) )
             CL = max(control.CL_range);
         end
    end
    
    CD = abs(CL) / control.CLCD;
    
    
    % generate output
    out.CL = CL;
    out.CD = CD;
    
end