function [ out ] = full_orbit(R0,V0,G,M,R_m,h_atm)
%Calculates the full orbit for selected initial conditions until sepcified
%end time
%   Might break off the orbit due to crash/flyby/end of time
    
    %%Input
    
    
    %%Functions
    [out_hk] = hyperbolic_kepler(R0,V0,G,M,R_m,h_atm);
    %while orbit
    %    [out_a] = in_atmosphere(??);
    %    [out_ek] = eliptic_kepler(??);
    %end
    %%Output
    out.R = 0;
    out.V = 0;
    out.V_rel = 0;
    out.a = 0;
    out.a_g = 0;
    out.a_aero = 0;
end

