function [ out ] = hyperbolic_kepler(R0,V0,A0,G,M,R_m,h_atm,dt_kep_init)
%Calculates the orbit from initial conitions until first entry into the atmosphere

%%Input


%%Functions
[out_param] = hk_Parameters(R0,V0,A0,G,M,dt_kep_init);
[out_end] = hk_Endpoint(out_param,G,M,R_m,h_atm);
%%Output
out = out_end;

end

