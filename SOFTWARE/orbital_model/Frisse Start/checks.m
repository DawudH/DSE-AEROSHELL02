function [ out_c ] = checks( R, V, t, tend, R_m, h_atm, G, M_mars, inatmos, crash_margin, round )
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
    
    %check if it makes at least one orbit.
    if round>=1
        out_c.orbit = true;
    else
        out_c.orbit = false;
    end
    
    % If the spacecraft is within the atmosphere boundary
    % or in orbit (we only compute the in atmos part then..) and not in
    % orbit without entering the atmosphere.. 
    if ( r <= (R_m + h_atm+100) ) || ( out_c.orbit && isreal(R) )
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
        else
            out_c.flyby = false;
        end
    else
        out_c.flyby = false;
    end
    
    


end

