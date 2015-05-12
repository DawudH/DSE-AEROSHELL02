clear
%DSE-02-CIA
%Created by: B.R. van Dongen
%Created on: 06-05-2015

%Description:
%
%
%
%
%
%

%--------------------------Orbit parameters--------------------------------
%--------------------------------------------------------------------------
Cd   = [1.2 1.5 1.5]                                     ;         %-
qmax = [1200 1000 200]                         ;        %Pa

mf = [2 2 1]                                    ; %mf = 1.0 for supersonic 
                                                  % and 2.0 for hypersonic

pstat = [100 100 100];
deltapstaticbar = pstat./qmax./Cd                    ;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%------------------------Gas characteristics-------------------------------
%--------------------------------------------------------------------------
M = 14                                          ;   %kg/kmole
R = 8314.472                                    ;   %J/kmole*K
T = 273.15                                      ;   %Gas temperature (in K)

molweight   = M*ones(1,3);
inflgastemp = T*ones(1,3);
infl_sys_mass_frac = 0.3;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%------------------------Geometric parameters------------------------------
%--------------------------------------------------------------------------
Do              = [3.0 15.0 50.0]               ;
Di              = 0.47*ones(1,3)                  ;
theta           = 60.0*ones(1,3)*pi./180.        ;
L               = 1                              ;   %Characteristic length

N               = [7 1 1]                        ;   %Number of toroids:

I               = [18 18 18]                     ;   %Number of radial straps
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%----------------------------Contingencies---------------------------------
%--------------------------------------------------------------------------
eta_p           = 1.25                           ; 
eta_g           = 1.25                           ;
eta_fiber       = 4.                             ;
eta_film        = 4.                             ;
eta_axial       = 4.                             ;
eta_toroid      = 4.                             ;
eta_radial      = 4.                             ;
eta_gores       = 4.                             ;

seam_allowance  = 1.05                           ;
gasb_allowance  = 1.05                           ;

mass_margin = 0.05                               ;
coating_mass_fraction = 0.5;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%------------------------------Materials-----------------------------------
%--------------------------------------------------------------------------
selected_material_gore = 'Kevlar 49'; %Gore material
selected_material_rad  = 'Kevlar 49'; %Radial strap material
selected_material_axi  = 'Kevlar 49'; %Axial strap material
selected_material_gasb = 'Kapton'; %Gas barrier material
selected_material_toro = 'Kevlar 49'; %Fiber or film material

TensYield_Nondim = [16583 ; 207000 ; 212400 ; 237453 ; 45000 ; 383918 ; 350999 ; 2220007 ; 36059 ; 231346; 22653];
Density = [1420;1440;1440;1700;1380;1540;970;1390;1470;1410;2700];
Strain_break = 1/100.*[72;3.6;2.4;1.4;30.5;3.5;3;4.4;42;4.3;0.];
Materials = { {'Kapton'}, {'Kevlar 29'}, {'Kevlar 49'}, {'M5'}, {'Nomex'}, {'PBO Zylon'}, {'Spectra 2000'}, {'Technora'}, {'Upilex-25S'},{'Vectran'},{'Aluminium'}};


%Toroid (fibres/film)
mat_idx_toro = find(strcmp([Materials{:}], selected_material_toro));
sigbar_toro = TensYield_Nondim(mat_idx_toro);
rho_toro    = Density(mat_idx_toro);
tmin_toro   = 2.54e-5;

typeoftoroid = 0; % 0 if Fibre ; 1 if Coat; 2 if Film

    %Fibres
beta            = 75/180*pi;
fiber_gap_ratio = 0.05;

mat_idx_gasb = find(strcmp([Materials{:}], selected_material_gasb));
sigbar_gasb = TensYield_Nondim(mat_idx_gasb);
rho_gasb = Density(mat_idx_gasb);
tmin_gasb = 5.08e-5;
 
%Axial
mat_idx_axi = find(strcmp([Materials{:}], selected_material_axi));
sigbar_axi = TensYield_Nondim(mat_idx_axi);
rho_axi    = Density(mat_idx_axi);    

%Radial
mat_idx_rad = find(strcmp([Materials{:}], selected_material_rad));
sigbar_rad = TensYield_Nondim(mat_idx_rad);
rho_rad    = Density(mat_idx_rad);    
    
%Gores
mat_idx_gore = find(strcmp([Materials{:}], selected_material_gore));
max_strain_gore = 1/100.*(10);
sigbar_gore = TensYield_Nondim(mat_idx_gore);
rho_gore        = Density(mat_idx_gore);
tmin_gore       = 1.4732e-5;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Define diameter ratio's
ksi_i       = Di./Do;
ksi_t       = [(1-ksi_i(1))./((2*N(1)-1)*sin(theta(1))+1-cos(theta(1))) 1/N(1) 1/N(1)];
ksi_d       = [ksi_i(1)+2*ksi_t(1) 1 1];


%Calculate toroid diameter
Dt          = ksi_t.*Do;

%Non-dimensionalize outer diameter
Dobar       = Do./L;

%Calculate Area Ratio (AR)
AR = [ 1-ksi_i(1).^2         1-ksi_i(2).^2       4*ksi_t(3).*(1-ksi_t(3)) ];

%Calculate dimensionless parameters
Cbar = N.*(1-ksi_t-ksi_t.*(N-1).*sin(theta));
Sbar = 4*pi*ksi_t./AR.*Cbar;
Vbar = Sbar./4;

%Calculate angles: 
%theta_c : undeformed half-cone constructed angle
%theta_t : torus attachment angle
%theta_d : deflected angle
%theta_h : heat-shield attachment angle
dzeta_all = AR./(4*mf.*ksi_t.*(1-ksi_t));
for i=1:3
    dzeta = dzeta_all(i); 
    thetadt0 = [60.65 60]*pi/180;
    options = optimoptions('fsolve','Display','off');
    F = @(thetadt) [tan(thetadt(1))+1/tan(2*thetadt(2))-1/sin(2*thetadt(2))*dzeta,
     sin(theta(i))/sin(thetadt(2))*(thetadt(1)-thetadt(2))-sin(thetadt(1)-thetadt(2))];
    [thetadt] = fsolve(F,thetadt0,options);
    
    theta_c(i)  = theta(i)               ;
    theta_d(i)  = thetadt(2)             ;
    theta_t(i)  = thetadt(1)             ;
    theta_h(i)  = 2*thetadt(2)-thetadt(1);
end
theta_h(1) = theta_c(1);
theta_d(1) = theta_c(1);
theta_t(1) = theta_c(1);


%Calculate minimum dimensionless gas pressure
pminbar = eta_p.*[AR(1).*tan(theta_c(1)).*sin(theta_c(1))/3./ksi_t(1) AR(2).*sin(theta_t(2))/pi./ksi_t(2).^2./cos(theta_h(2)) AR(3).*sin(theta_t(3))/pi./ksi_t(3).^2./cos(theta_h(3))];

%Calculate inflation gas mass
pbar = pminbar + deltapstaticbar;
Gbar = 9.80665.*L./(R./M.*T); 
mgasbar = Gbar.*pbar.*Vbar.*ksi_t.*Dobar.*eta_g;

minflsysbar = mgasbar*infl_sys_mass_frac;

%Area of IAD:
A_IAD = [pi/4*(Do(1).^2-Di(1).^2) pi/4*(Do(2).^2-Di(2).^2) pi/4*(Do(3).^2-(Do(3)-2*Dt(3)).^2)];
%or:
A_IAD = AR.*pi/4.*Do.^2;
    
mfactor = A_IAD.*Cd.*qmax./9.80665;

%--------------------------Reinforcing elements----------------------------
%--------------------------------------------------------------------------

%Radial strap length:
Lrbar = [ksi_t(1).*(2*(N(1)-1)+pi) (1-ksi_i(2)-ksi_t(2)+cos(theta(2))*ksi_t(2))/(2*sin(theta(2)))+pi.*ksi_t(2) (1-ksi_i(3)-ksi_t(3))/(2*sin(theta(3)))+ksi_t(3)*(3*pi./4+1)];

%Gores
Agorebar = 1./sin(theta_c)+2.*pi.*ksi_t.*(1-ksi_t)./AR;

%Gores max curvature
rcbar = pi./I*(4*max_strain_gore.^2+1)./(8*max_strain_gore);
rc    = rcbar.*Do;

%Radial strap mass
mradialbar = 1./sigbar_rad.*Lrbar./cos(theta_h).*Dobar.*eta_radial;

S = Sbar.*A_IAD;
a = rho_toro.*S.*tmin_toro./mfactor;
b = 1./sigbar_toro .* pminbar .* Sbar .* Dobar .* eta_toroid .* 1./4 .* ksi_t .* (2 - 3.*ksi_t) ./ (ksi_d - 2 .* ksi_t);

if typeoftoroid==0 %Fibre reinforced
    
    %Component 1 -> Fiber:
    mfilmbar       = max(a,1./sigbar_toro.*(1+1./(tan(beta).^2)).*pminbar.*Sbar.*ksi_t./2.*Dobar.*eta_fiber);

    %Component 2 -> Gas barrier:
    mgasbarrierbar  = max(rho_gasb.*S.*tmin_gasb./mfactor,b*fiber_gap_ratio);
   
    %Component 3 -> Axial straps:
    maxialbar = 1./sigbar_axi.*pminbar.*Sbar.*ksi_t./4.*Dobar.*eta_axial;
    
    %Component 4 -> Coating:
    mcoatingbar = 0.5*mfilmbar;

    
    mtoroidbar = mfilmbar+mgasbarrierbar+maxialbar+mcoatingbar;
   
elseif typeoftoroid==1 %Coated
    
    %Component 1 -> Film:
    mfilmbar = max(a,b);
    
    %Artificial adjustment
    mfilmbar(1)=mfilmbar(1)./4;
   
    %Component 2 -> Coating:
    mcoatingbar = 0.5*mfilmbar;
    
    maxialbar = [nan nan nan];
    mgasbarrierbar = [nan nan nan];
    
    mtoroidbar = mfilmbar + mcoatingbar;

elseif typeoftoroid==2 %Film
    mfilmbar = max(a,b);
    
    %Artificial adjustment
    mfilmbar(1)=mfilmbar(1)./4;
   
    maxialbar = [nan nan nan];
    mgasbarrierbar = [nan nan nan];
    mcoatingbar = [nan nan nan];
    mtoroidbar = mfilmbar;
end

%Dimensionless gore mass
mgoresbar = max(rho_gore.*S.*tmin_gore./mfactor, ...
                1./sigbar_gore.*rcbar.*Agorebar.*Dobar.*eta_gores);
mgoresbar(3) = 0;

mgoresbar = seam_allowance*mgoresbar;
mgasbarrierbar = gasb_allowance*mgasbarrierbar;

mmarginbar = mass_margin*ones(1,3);
msumbar = mgasbar+minflsysbar+mtoroidbar+mradialbar+mgoresbar;
mtotalbar = msumbar./(1-mmarginbar);
mtotal = mtotalbar.*mfactor;
pmin = pminbar.*qmax.*Cd;
%Tabulate results:
wline = [nan nan nan];
data = [mf;Cd;pstat;qmax;wline;N;I;Do;Di;Dt;ksi_t;AR;Vbar;Sbar;Cbar;theta_c*180./pi;theta_d*180./pi;theta_t*180./pi;theta_h*180./pi;
     Lrbar;Agorebar;rcbar;wline;[1 1 1];[2 2 2]; [3 3 3]; [4 4 4];
     wline;molweight;inflgastemp;pminbar;pmin;mgasbar;minflsysbar;wline;mfilmbar;mcoatingbar;mgasbarrierbar;maxialbar;mradialbar;mgoresbar;wline;
     mmarginbar;mtotalbar;wline;mfactor;mtotal];
colheadings = {'Stacked toroid', 'Tension cone', 'Trailing IAD'};
rowheadings = {'Regime','Drag coefficient [-]','Static pressure [Pa]','Dynamic pressure [Pa]',' ','Number of toroids [-]','Number of radial straps [-]', 'Outer diameter [m]', 'Inner diameter [m]', 'Toroid diameter [m]','Ratio of minor toroid to outer diameter [-]',... 
              'Projected area ratio [-]', 'Toroid volume [-]', 'Toroid surface area [-]', 'Toroid circumference [-]',...
              'Constructed angle [deg]', 'Deflected angle [deg]', 'Torus attachment angle [deg]', 'Heat shield attachment angle [deg]',...
              'Radial strap length [-]','Gores surface area [-]','Gores max curvature [-]',' ',...
              'Fibre/Film material', 'Gore material','Axial strap material','Radial strap material', ' ',...
              'Inflation gas molecular weight [kg/mole]','Inflation gas temperature [K]','Inflation pressure [-]','Inflation pressure [Pa]', 'Inflation gas mass [-]','Inflation system mass [-]',' ',...
              'Toroid fiber/film mass [-]','Toroid coating mass [-]','Toroid gas barrier mass [-]', 'Axial straps mass [-]', 'Radial straps mass [-]', 'Gores mass [-]',' ',...
              'Mass margin [-]','Total mass [-]',' ','Mass Factor [kg]', 'Total mass [kg]'};        
wid = 16;
%fms = {'0.5E','0.5E','0.5E'};
fms = {'d','d','d'};
displaytable(data, colheadings, wid,fms,rowheadings)

