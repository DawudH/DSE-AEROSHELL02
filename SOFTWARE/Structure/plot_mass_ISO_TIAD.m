clear all
close all
clc


addpath('.\matlab2tikz')

figure(1)

Do=12
q=10000
model=1
Material='Vectran'
Cd=1.04

plotnr= 2

run Mass_est_Isotensoid_function.m

ml=[];
ql=[];
Dl=[];
Al=[];
BCl=[];

mm=[];
qm=[];
Dm=[];
Am=[];
BCm=[];

if model==1
    mname='\textbf{Isotensoid}';
end
if model==2
    mname='\textbf{Trailing IAD}';
end
if plotnr==1
    for Do=5:1:15
        run Mass_est_Isotensoid_function.m
        ml=[ml; m]
        ql=[ql; q]
        Dl=[Dl; Do]
        Al=[Al; A]
        BCl=[BCl; BC] 
    end
        plot(Dl,BCl)
        title({mname},'Interpreter','LaTex','FontSize',13)
        legend(Material,'Interpreter','LaTex','FontSize',12,2)
        xlabel('Dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
        ylabel('Mass[kg]','Interpreter','LaTex','FontSize',14)
        grid on
end

if plotnr==2
    for q=1000:500:10000
       
        
        mm=[mm ml]
        qm=[qm ql];
        Dm=[Dm Dl];
        Am=[Am Al];
        BCm=[BCm BCl];
        
        
        ml=[]
        ql=[]
        Dl=[]
        Al=[]
        BCl=[]

        for Do=5:0.5:15
        run Mass_est_Isotensoid_function.m
        ml=[ml; m];
        ql=[ql; q];
        Dl=[Do; Dl];
        Al=[Al; A];
        BCl=[BCl; BC];
        end  
    end
        y=qm(1,:)'
        x=Dm(end:-1:1,1)
        z=mm'
        contourf(x,y,z)
        colormap parula
        h = colorbar
        title({mname},'Interpreter','LaTex','FontSize',13)
        ylabel('Dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
        xlabel('Diameter [m]','Interpreter','LaTex','FontSize',14)
        ylabel(h,'Decelerator Mass [kg]','Interpreter','LaTex','FontSize',14)
        
        grid on
        yl = ylim;
        xl = xlim;
        %matlab2tikz('.\LaTeX\test.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
        %fix_grid_contourf('.\LaTeX\test.tikz')
end



if plotnr==3
    markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
    Materials_list={'Vectran','Kevlar 29','Kevlar 49','Upilex-25S', 'Nomex', 'Nextel 610'}
    xrange=size(Materials_list)
    xrange=xrange(1,2)
    for i= 1:xrange;
        Material=Materials_list{i}
        ml=[];
        ql=[];
        Dl=[];
        Al=[];
        BCl=[];
        
        for q=500:500:5000
            run Mass_est_Isotensoid_function.m
            ml=[ml; m];
            ql=[ql; q];
            Dl=[Do; Dl];
            Al=[Al; A];
            BCl=[BCl; BC]; 
    end
        plot(ql,BCl,markers{i+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    end
        
        title({mname},'Interpreter','LaTex','FontSize',13)
        legend(Materials_list,'Interpreter','LaTex','FontSize',12,2)
        xlabel('Peak dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
        ylabel('Ballistic coefficient[$\frac{kg}{m^2}$]','Interpreter','LaTex','FontSize',14)
        grid on
        %matlab2tikz('.\LaTeX\Tikz\ISO_mat.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
        
end       
        

    
    
    
    
    
    