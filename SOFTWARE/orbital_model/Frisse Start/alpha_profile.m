function out = alpha_profile(t,aero_coef)

    Alpha = [20]*pi/180;
    t_change = [226];
    
    t_check = t_change >= t;
    if sum(t_check) == 0
        alpha = Alpha(end);
    else
        alpha = Alpha(find(t_check,1));
    end

    [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);

    out.CLA = CLA;
    out.CDA = CDA;
    out.CMYA = CMYA;
    out.alpha = alpha;
    out.error = 0;
    out.error_I = 0;
end