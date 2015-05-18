
%User defined vatiables
Do =Do;
A= pi/4*Do^2;
q=q;

model=model;    
Material=Material;

%dfml=[0.0916;0.1457;0.2080;0.1810;0.3750;0.3778;0.4001;0.2780];
dfml=[          0.1457;          0.1457;         0.2080;     0.1810;         0.3750;             0.3778;     0.4001;         0.2780];
sigy=[           3.2;                3.2;           2.92;     3.00;              3.00;            0.52;        0.61;           3.1;           ]
rho	=[          1410;           1410;              1440;        1440;          1440;                1470;        1410;           3900;             ]        
Materials = { {'Vectran'}, {'Coated Vectran'}, {'Kevlar 29'}, {'Kevlar 49'}, {'Coated Kevlar'},{'Upilex-25S'}, {'Nomex'}, {'Nextel 610'},};
mat_idx = find(strcmp([Materials{:}], Material));
kfl=(sigy*10^9)./rho
dfm=dfml(mat_idx);


qbar_q  =1;
kc      =299000;
%kf      =98500;
kf=kfl(mat_idx)
KD      =2.9;
KC      =1.25;



if model==1
    %Given variables ISO tensoid
    %Cd      =1.06;
    Cd=Cd
    P_q     =2;
    Tbar    =0.44;
    Af_piR  =2.81;
    fbar    =0.09;
    eta     =0.05;
    lm_R    =2.56;
    lt_R    =0;
    xi      =0.40;
    beta    =2;

    b       =1.1e-5;
    c       =4.02;
end
if model==2
    %Given variables Balute
    %Cd      =0.63;
    Cd=Cd
    P_q     =2.4;
    Tbar    =0.865;
    Af_piR  =4;
    fbar    =0.135;
    eta     =0.1;
    lm_R    =pi;
    lt_R    =4;
    xi      =0;
    beta    =2;

    b       =6.9e-5;
    c       =7.41;
end


%Mass computation
qCdA05=q*sqrt(Cd*A);
df= max((KD*beta*fbar*P_q/(2*(1+2*eta)*kf*sqrt(Cd*pi))) * qCdA05 , dfm)*5;
BC= (KD/ (sqrt(Cd*pi)*kc*(1-xi^2)^1.5)) *  (P_q/Cd*lm_R*Tbar/((1+2*eta)^3) + lt_R)  * (qCdA05)  + (KC/(Cd*(1+2*eta)^2)) * (Af_piR + 4*pi*eta*(1+eta)) * (df/(1-xi^2))
%BC= b*qCdA05 + c*df
m=Cd*A*BC

%Imperial coefficients
BCimp=BC*0.0929/0.4536;













