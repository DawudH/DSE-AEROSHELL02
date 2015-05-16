%User defined vatiables
Do =12;
A= pi/4*Do^2;
q=3000;

model=1
Material='Vectran'

dfml=[0.0916;0.1457;0.2080;0.1810;0.3750;0.3778;0.4001;0.2780];
Materials = { {'Vectran'}, {'Coated Vectran'}, {'Kevlar 29'}, {'Kevlar 49'}, {'Coated Kevlar'},{'Upilex-25S'}, {'Nomex'}, {'Nextel 610'}};
mat_idx = find(strcmp([Materials{:}], Material));

dfm=dfml(mat_idx)

qbar_q  =1;
kc      =299000;
kf      =98500;
KD      =2.9;
KC      =1.25;

BCl=[];
BCimpl=[];
ql=[];
qimpl=[];




for q=20:1:4000
    
for i=[1 2]
    if i==1
        %Given variables ISO tensoid
        Cd      =1.06;
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
    if i==2
        %Given variables Balute
        Cd      =0.63;
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

    qCdA05=q*sqrt(Cd*A)

    %Mass computation
    df= max((KD*beta*fbar*P_q/(2*(1+2*eta)*kf*sqrt(Cd*pi))) * qCdA05 , dfm)
    BC(i)= (KD/ (sqrt(Cd*pi)*kc*(1-xi^2)^1.5)) *  (P_q/Cd*lm_R*Tbar/((1+2*eta)^3) + lt_R)  * (qCdA05)  + (KC/(Cd*(1+2*eta)^2)) * (Af_piR + 4*pi*eta*(1+eta)) * (df/(1-xi^2))
    %BC(i)= b*qCdA05 + c*df
    m(i)=Cd*A*BC(i)

    %Imperial coefficients
    BCimp(i)=BC(i)*0.0929/0.4536;

end

data = [BC; BCimp; m];
colheadings = {'Isotensoid', 'Ballute'};
rowheadings = {'BC(SI)' 'BC(Imperial)', 'mass'};     
wid = 16;
fms = {'d','d'};
displaytable(data, colheadings, wid,fms,rowheadings)


BCl=[BCl; BC];
BCimpl=[BCimpl; BCimp];
ql=[ql; qCdA05];
qimpl=[qimpl; qCdA05*0.3048/4.448];
end



loglog(qimpl,BCimpl(:,model))









