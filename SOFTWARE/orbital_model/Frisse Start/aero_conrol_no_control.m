function [out] = aero_conrol_no_control(state,control,aero_coef,V,R,q,Omega_m,m,g_earth)

        if q > 0
        
%             % determine range of alpha
%             alpha_range = min(control.alpha_range) : 0.1*pi/180 : max(control.alpha_range);
%             [CLA_range, CDA_range, CMYA] = aero_coef.aeroCoeffs(alpha_range);
% 
%             %calculate the velocity of the atmosphere
%             Vatm = cross(Omega_m,R);
%             %calculate the direction of the velocity of the s/c wrt the atmosphere
%             vel_unit = (V - Vatm) / norm(V - Vatm);
% 
% 
%             a_human_range_mag = zeros(length(CDA_range),1);
%             for i = 1:length(CDA_range)
%                 a_human_range_mag(i) = norm(cross(vel_unit,[0,0,1]) * CLA_range(i) * q / m - vel_unit * CDA_range(i) * q / m);
%             end
% 
%             % find value closest to the setpoint
%             errors = a_human_range_mag - control.a;
%             if errors > 0 % all positive, so control.a is smaller then all possible setpoints
%                 % maximize lift to get out of the atmosphere..
%                 % max lift at min alpha...
%                 alpha = min(alpha_range);
%             elseif errors < 0 % all negative, so control.a is smaller then all possible setpoints
%                 % minimize lift to get into atmosphere quick..
%                 % min lift at max alpha
%                 alpha = max(alpha_range);
%             else % combination of positive and negative values, so an optimal point can be found!
%                 % problem is that the lift should be negative when below
%                 % the setpoint accel and positive when above setpoint accel
%                 % so if above the setpoint lift positive -> alpha
%                 % negative..
%                 
%                 if state.a > control.a
%                     % negative alpha
%                     index = find(alpha_range < 0);
%                 else
%                     % positive alpha
%                     index = find(alpha_range > 0);
%                 end
%                 
%                 [~, index_a] = min(abs(errors(index)));
%                 loc = index(index_a);
%                 %disp(['closest to setpoint possible a_human: ' num2str(a_human_range_mag(loc) ) ])
%                 alpha = alpha_range(loc);
%             end
%             
% 
%             
%         else
%             alpha = state.alpha;


            if state.a < control.a
            % If state.a < control.a, decrease CL to allow higer acceleration
                 alpha = state.alpha - control.dalpha;

            elseif state.a == control.a
            % If state.a = control.a, keep the same CL
                alpha = state.alpha;

            else
            % If state.a > control.a, increase CL to allow lower acceleration
                alpha = state.alpha + 3*control.dalpha;
                 
            end

            % check if alpha is out of bounds
                 if ( alpha > max(control.alpha_range) )
                     alpha = max(control.alpha_range);
                 elseif ( alpha < min(control.alpha_range) )
                     alpha = min(control.alpha_range);
                 end


        else
            alpha = state.alpha;

        end
        
        [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);

    % BEUN
    e = 0;
    error_I = 0;
            
    % generate output
    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    out.error = e;
    out.error_I = error_I;
end