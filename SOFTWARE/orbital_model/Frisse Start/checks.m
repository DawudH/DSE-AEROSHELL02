function [ out_c ] = checks( R, V, t, tend, R_m, h_atm, G, M_mars, inatmos, crash_margin )
%CHECKS Summary of this function goes here
%   Detailed explanation goes here

    % calculate the magnitude of the distance
    r = norm(R);
    
    % If mission is at the end of time ;)
    if ( t >= tend )
        out_c.t_end = true;
    else
        out_c.t_end = false;
    end

    % If the spacecraft is within the atmosphere boundary
    if ( r < (R_m + h_atm) ) 
        out_c.in_atmos = true;
    else
        out_c.in_atmos = false;
    end
    
    % if the spacecraft is crashing on the surface + some margin
    if ( r < (R_m + crash_margin) )
        out_c.crash = true;
    else
        out_c.crash = false;
    end
    
    % if R is imaginary (hyperbola does noet reach atmosphere)
    % there is a flyby!
    % the other condition is when the Vesc is higher then the velocity at
    % the end of the atmosphere
    % escape velocity at the inital point
     
    if (isreal(R) == false)
        out_c.flyby = true;
    elseif (inatmos && (out_c.in_atmos == false))
        
        V_esc = sqrt(G*M_mars * 2 / r);
        
        if (norm(V) > V_esc)
            out_c.flyby = true;
        end
    else
        out_c.flyby = false;
    end
    
    


end

