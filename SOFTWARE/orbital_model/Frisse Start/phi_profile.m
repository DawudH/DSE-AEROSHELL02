function out = phi_profile(t)


%     phi = [0 30 -30 30 -30 0 20 -20 20]*pi/180;
%     t_change = [200 270 350 450 550 1200 1300 1400];
    
    phi = [0 35]*pi/180;
    t_change = [160];

    
    t_check = t_change >= t;
    if sum(t_check) == 0
        check = length(phi);
    else
        check = find(t_check,1);
    end
    
    out = phi(check);
end