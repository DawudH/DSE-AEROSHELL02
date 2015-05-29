function out = alpha_profile(t,aero_coef,control,state)

    Alpha = [10 17 16 15 10]*pi/180;
    t_change = [226 234 250 270];
    
    t_check = t_change >= t;
    if sum(t_check) == 0
        check = length(Alpha);
    else
        check = find(t_check,1);
    end
    
    if (state.alpha - control.dalpha) < Alpha(check)
        alpha = state.alpha - control.dalpha;
    elseif (state.alpha + control.dalpha) > Alpha(check)
        alpha = state.alpha + control.dalpha;
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