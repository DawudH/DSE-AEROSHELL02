function out = phi_profile(t)


%     phi = [0 30 -30 30 -30 0 20 -20 20]*pi/180;
%     t_change = [200 270 350 450 550 1200 1300 1400];
%     %d=12 one-time
%     phi = [0 32 0 32]*pi/180;
%     t_change = [160 250 350];
%     
%     %d=6 one-time
%     phi = [0 32 0 32]*pi/180;
%     t_change = [160 270 400];
%     

	% Aero capture
    phi = [0 38.471 0 0]*pi/180;
    t_change = [240 1000 1100];   
            
            % density * 1.1
%             phi = [0 52.73 0 0]*pi/180;
%             t_change = [240 1000 1100];   

            % density * 0.90
%             phi = [0 7.2 0 0]*pi/180;
%             t_change = [240 1000 1100];   
    
%     % second orbit
%     phi = [129 60 30]*pi/180;
%     t_change = [500 700]; 
    
            % density * 1.1
%             phi = [130 80 30]*pi/180;
%             t_change = [500 700]; 

            % density * 0.90
%             phi = [129 0 30]*pi/180;
%             t_change = [500 700];   

    t_check = t_change >= t;
    if sum(t_check) == 0
        check = length(phi);
    else
        check = find(t_check,1);
    end
    
    out = phi(check);
end