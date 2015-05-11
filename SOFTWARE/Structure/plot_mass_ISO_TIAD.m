addpath('.\matlab2tikz')

Do=12
q=3000
model=1
Material='Vectran'

plotnr= 2

run Mass_est_Isotensoid_function.m

ml=[]
ql=[]
Dl=[]
Al=[]
BCl=[]

mm=[]
qm=[]
Dm=[]
Am=[]
BCm=[]

if model==1
    mname='\textbf{Isotensoid}'
end
if model==2
    mname='\textbf{Trailing IAD}'
end
if plotnr==1
    for q=500:500:5000
        run Mass_est_Isotensoid_function.m
        ml=[ml; m]
        ql=[ql; q]
        Dl=[Do; Dl]
        Al=[Al; A]
        BCl=[BCl; BC] 
    end
        plot(ql,ml)
        title({mname},'Interpreter','LaTex','FontSize',13)
        legend(Material,'Interpreter','LaTex','FontSize',12,2)
        xlabel('Dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
        ylabel('Mass[kg]','Interpreter','LaTex','FontSize',14)
        grid on
end

if plotnr==2
    for q=500:500:5000
       
        
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

        for Do=5:1:12
        run Mass_est_Isotensoid_function.m
        ml=[ml; m]
        ql=[ql; q];
        Dl=[Do; Dl];
        Al=[Al; A];
        BCl=[BCl; BC];
        end  
    end
        x=qm(1,:)'
        y=Dm(end:-1:1,1)
        z=mm
        figure(1)
        contourf(x,y,z)
        colormap parula
        h = colorbar
        title({mname},'Interpreter','LaTex','FontSize',13)
        xlabel('Dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
        ylabel('Diameter [m]','Interpreter','LaTex','FontSize',14)
        ylabel(h,'Decelerator Mass [kg]','Interpreter','LaTex','FontSize',14)
        
        grid on
        yl = ylim;
        xl = xlim;
        matlab2tikz('.\LaTeX\test.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
        fix_grid_contourf('.\LaTeX\test.tikz')

        
end
    
    
    
    
    
    
    