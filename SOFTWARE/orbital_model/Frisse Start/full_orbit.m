function [ out ] = full_orbit(R0, V0, A0, G, M_mars, R_m, h_atm, atm, dt_kep_init, dt_atmos, m, Omega_m, S, control, tend, crash_margin, g_earth, aero_coef, use_control, multiple_orbits,use_alpha_profile,alpha_init,dalphadt)
%Calculates the full orbit for selected initial conditions until sepcified
%end time
switch nargin
    case 22
        control.alpha_init = alpha_init;
        no_waitbar = true;
    case 23
        control.alpha_init = alpha_init;
        control.dalphadt = dalphadt;
        control.Kp = control.Kp/control.dalpha*control.dalphadt*dt_atmos; % proportional gain
        control.Ki = control.Ki/control.dalpha*control.dalphadt*dt_atmos; % integral gain
        control.Kd = control.Kd/control.dalpha*control.dalphadt*dt_atmos; % differential gain
        control.dalpha = control.dalphadt*dt_atmos;
        
        no_waitbar = true;
    otherwise
        no_waitbar = true;
end
    
if no_waitbar == false
    h = waitbar(0,'Initializing waitbar...');
end
    %Condition to run the script
    orbit = true;

    %Counters
    tp(1) = 0; %Time for plot
    t = 0; %Time
    round = 0; %Number of orbits around mars
    i = 1; %Number of loops
    
    % give initial control state
    alpha = control.alpha_init;
    [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
    
    %%Functions
    [out_hk] = hyperbolic_kepler(R0,V0,A0,G,M_mars,R_m,h_atm,dt_kep_init);
    %%Inputs for while loop
    R(1,:) = out_hk.end.R;
    speed_sound(1,:) = out_hk.end.speed_sound;
    V(1,:) = out_hk.end.V;
    M(1,:) = out_hk.end.M;
    speed_sound(1,:) = out_hk.end.speed_sound;
    A(1,:) = out_hk.end.A;
    Ag(1,:) = out_hk.end.Ag;
    Ad(1,:) = out_hk.end.Ad;
    Al(1,:) = out_hk.end.Al;
    J(1,:) = out_hk.end.J;
    q(1,:) = out_hk.end.q;
    T(1,:) = out_hk.end.T;
    rho(1,:) = out_hk.end.rho;
    Alpha(1) = alpha;
    dAlpha_dt(1) = 0;
    CD(1) = CDA / S;
    CL(1) = CLA / S;
    
    a_prev = A(1,:);
    
    state.CL = CL(1);
    state.CD = CD(1);
    state.a = norm(A(1,:) - Ag(1,:));
    state.h = norm(R(1,:)) - R_m;
    state.alpha = alpha;
            
    %Get initial values for conditions
    [out_c] = checks( R(1,:), V(1,:), t, tend, R_m, h_atm, G, M_mars,g_earth, false, crash_margin, false, round, state, control );
    out_c.in_atmos = true;
    
    

    %%Functions
    %As long as the s/c is in orbit keep calculating the next position
    while orbit
        %When the s/c is in the atmosphere use the numerical solution including aerodynamic forces
        if out_c.in_atmos
             % determine new cl and cd param
                 % start controlling once the accel is above 1.5g
                 if use_control && out_c.use_control
                         [aero_param] = aero_conrol(state,control,aero_coef,dt_atmos);
                         CL(i+1) = aero_param.CLA / S;
                         CD(i+1) = aero_param.CDA / S;
                         alpha = aero_param.alpha;
                         control.error = aero_param.error;
                         control.error_I = aero_param.error_I;
                 else
                        control.error = 0;
                        control.error_I = 0;
                        CL(i+1) = CL(i);
                        CD(i+1) = CD(i); 
                 end
                 
                 if use_alpha_profile && ~use_control
                     aero_param = alpha_profile(t,aero_coef,control,state);
                     CL(i+1) = aero_param.CLA / S;
                     CD(i+1) = aero_param.CDA / S;
                     alpha = aero_param.alpha;
                 end
                
            [out_o] = in_atmosphere( V(i,:), R(i,:), A(i,:), a_prev, J(i,:), atm, CL(i+1), CD(i+1), dt_atmos, R_m, Omega_m, S, m );
            %update state
            state.CL = CL(i);
            state.CD = CD(i);
            state.a = norm(A(i,:) - Ag(i,:));
            state.h = norm(R(i,:)) - R_m;
            state.alpha = alpha;
            a_prev = A(i,:);
            t = t + dt_atmos;
        %When the s/c is not in the atmosphere use a kepler orbit
        else
            orbit_init.R = R(i,:);
            orbit_init.V = V(i,:);
            orbit_init.a = A(i,:);
            [out_e] = eliptic_kepler(R(i,:),V(i,:),A(i,:),G,M_mars,dt_kep_init,orbit_init);
            out_o = out_e;
            round = round + 1;
            a_prev = out_o.A;
            tkep = tp(i);
            t = t + out_o.t_kep;
            %readjust attitude to get back to initial alpha
            alpha = control.alpha_init;
            [CLA, CDA, CMYA] = aero_coef.aeroCoeffs(alpha);
            CL(i+1) = CLA / S;
            CD(i+1) = CDA / S;
            state.alpha = alpha;
            control.error = 0;
            control.error_I = 0;

        end
        
        %Function to check when to end the orbit              
        [out_c] = checks( out_o.R, out_o.V, t, tend, R_m, h_atm, G, M_mars,g_earth, out_c.in_atmos, crash_margin, out_c.use_control,round, state, control );
            
        %%New Inputs for while loop
        tp(i+1) = tp(i) + dt_atmos;
        R(i+1,:) = out_o.R;
        speed_sound(i+1,:) = out_o.speed_sound;
        V(i+1,:) = out_o.V;
        speed_sound(i+1,:) = out_o.speed_sound;
        M(i+1,:) = out_o.M;
        A(i+1,:) = out_o.A;
        Ag(i+1,:) = out_o.Ag;
        Ad(i+1,:) = out_o.Ad;
        Al(i+1,:) = out_o.Al;
        J(i+1,:) = out_o.J;
        q(i+1,:) = out_o.q;
        T(i+1,:) = out_o.T;
        rho(i+1,:) = out_o.rho;
        Alpha(i+1) = state.alpha;
        if i < 1
               dAlpha_dt(i+1) = 0;
        else
               %dAlpha_dt(i+1) = (Alpha(i-1) - 4*Alpha(i) + 3*Alpha(i+1))/(dt_atmos * 2);
               dAlpha_dt(i+1) = (Alpha(i+1) - Alpha(i))/(dt_atmos);
        end
        
        if out_c.crash || out_c.flyby || out_c.t_end || ( (multiple_orbits == false) && out_c.orbit)
            orbit = false;
        end
        i = i+1;
        if no_waitbar == false
            waitbar(t/tend,h,sprintf('%5.1f %...',t/tend*100))
        end
    end
    
    
    
    %%Output
    out.R = R;
    out.V = V;
    out.A = A;
    out.Ag = Ag;
    out.Ad = Ad;
    out.Al = Al;
    out.J = J;
    out.q = q;
    out.tp = tp;
    out.M = M;
    out.t = t;
    out.theta_p = out_hk.param.theta_p;
    out.rp = out_hk.param.rp;
    out.ra = out_hk.param.ra;
    out.theta0 = out_hk.param.theta;
    out.theta = out_hk.end.theta;
    out.a = out_hk.param.a;
    out.e = out_hk.param.e;
    out.rc = out_hk.end.rc;
    out.c = out_c;
    out.speed_sound = speed_sound;
    out.T = T;
    out.rho = rho;
    out.alpha = Alpha;
        a_human_mag = sqrt((out.Ad(:,1)+out.Al(:,1)).^2 + (out.Ad(:,2)+out.Al(:,2)).^2 + (out.Ad(:,3)+out.Al(:,3)).^2);
        maxaccel = max(a_human_mag)/g_earth;
    out.a_human_mag = a_human_mag;
    out.maxaccel = maxaccel;
    out.CL = CL;
    out.CD = CD;
    out.Kp = control.Kp;
    out.Ki = control.Ki;
    out.Kd = control.Kd;
    out.dAlpha_dt = dAlpha_dt;
    
    if exist('tkep','var')
        out.tkep = tkep;
        out.okep = out_e.param;
    end
    
    % output text
    

    disp(['rx = ' num2str(R0(1)) ' [m], CD = ' num2str(CDA / S) ' [-], CL = ' num2str(CLA / S) ' [-], in atmosphere: ' num2str(out_c.in_atmos) ', crashed: ' num2str(out_c.crash) ', in orbit: ' num2str(out_c.orbit) ', flyby: ' num2str(out_c.flyby) ', acceleration: ' num2str(maxaccel) ', time pased: ' num2str(t/(3600*24)) ' days' ])
    
    if no_waitbar == false
        %close waitbar
        close(h);
    end
end

