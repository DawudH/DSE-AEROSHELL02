function [ out ] = eliptic_kepler(R0,V0,A0,G,M,dt_kep_init,orbit_init)
%Calculates the orbit while in eliptic orbit outside atmosphere

%%Input


%%Functions
[out_param] = ek_Parameters(R0,V0,A0,G,M,dt_kep_init);
[out_end] = ek_Endpoint(out_param,orbit_init);
%%Output
out = out_end;
out.param = out_param;
out.V_aero = 0;
end
