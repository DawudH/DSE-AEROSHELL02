function out = alpha_profile(t,aero_coef,control,state,dt)
    %something else
    Alpha = [18 5 10 10]*pi/180;
    Alpha = [20 20 20 20]*pi/180;
    t_change = [1000 1100 1300];
    dalpha = -[0.5 0.5 0.5 0.5]*pi/180*dt;
    
    %one time under 3g
    Alpha = [10 10 30 0]*pi/180;
    t_change = [1000 1290 1385];
    dalpha = -[0.5 0.5 0.5 0.5 0.5]*pi/180*dt;
    
    t_check = t_change >= t;
    if sum(t_check) == 0
        check = length(Alpha);
    else
        check = find(t_check,1);
    end
    
    if (state.alpha - control.dalpha) < Alpha(check)
        alpha = state.alpha - dalpha(check);
    elseif (state.alpha + control.dalpha) > Alpha(check)
        alpha = state.alpha + dalpha(check);
    else
        alpha = Alpha(check);
    end
    [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);

    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    out.error = 0;
    out.error_I = 0;
end