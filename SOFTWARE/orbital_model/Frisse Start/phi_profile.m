function out = phi_profile(t)


%     phi = [0 30 -30 30 -30 0 20 -20 20]*pi/180;
%     t_change = [200 270 350 450 550 1200 1300 1400];

	%something else
    phi = [0 35 0 0]*pi/180;
    t_change = [160 1000 1300];    
    %d=12 one-time
    phi = [0 32 10 32]*pi/180;
    t_change = [160 250 350];
    
    %d=6 one-time
    phi = [0 32 0 32]*pi/180;
    t_change = [160 270 400];
    
    t_check = t_change >= t;
    if sum(t_check) == 0
        check = length(phi);
    else
        check = find(t_check,1);
    end
    
    out = phi(check);
end