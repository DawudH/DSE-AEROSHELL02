function out = alpha_profile(t,aero_coef,control,state,dt)
    %something else
    Alpha = [18 5 10 10]*pi/180;
    Alpha = [20 20 20 20]*pi/180;
    t_change = [1000 1100 1300];
    dalpha = -[0.5 0.5 0.5 0.5]*pi/180*dt;
    
    %one time under 3g
    Alpha = [13.2 33.2 13.2]*pi/180;
    t_change = [500 520];
    dalpha = -[0.5 0.5 0.5 0.5 0.5]*pi/180*dt;
    
    %Alpha = [9 30 9 9]*pi/180;
    
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
    [CL, CD, CMY] = aero_coef.aeroCoeffs(alpha);

    out.CL = CL;
    out.CD = CD;
    out.CMY = CMY;
    out.alpha = alpha;
    out.error = 0;
    out.error_I = 0;
end