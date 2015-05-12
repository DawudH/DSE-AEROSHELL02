function [ out ] = eliptic_kepler(R0,V0,A0,G,M,R_m,h_atm,dt_kep_init)
%Calculates the orbit while in eliptic orbit outside atmosphere

%%Input


%%Functions
[out_param] = ek_Parameters(R0,V0,A0,G,M,dt_kep_init);
[out_end] = ek_Endpoint(out_param,G,M,R_m,h_atm);
%%Output
out = out_end;

end
