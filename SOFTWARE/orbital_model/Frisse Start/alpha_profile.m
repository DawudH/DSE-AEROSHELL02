function out = alpha_profile(t,aero_coef,control,state)

    Alpha = [20 19 18 17 18 19 20 21 20 19 18 17]*pi/180;
    t_change = [245 250 255 260 265 270 275 280 290 300 310];
    
    Alpha = [20 5]*pi/180;
    t_change = [210];
    
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