m0=10000


mstructl=[]
mheatl=[]
Ql=[]
ql=[]


for q=1000:100:10000


mstruct=(0.0232*q ^0.1708)*m0;
mstructl=[mstructl mstruct];
ql=[ql q];


end

figure(1)
plot(ql,mstructl,'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
xlabel('Dynamic pressure [Pa]','Interpreter','LaTex','FontSize',14)
ylabel('Forebody sttructural mass[kg]','Interpreter','LaTex','FontSize',14)
grid on



for Q=1000:1000:100000
    
mheat=(0.00091*Q^0.51575)*m0;
mheatl=[mheatl mheat];
Ql=[Ql Q];

end

figure(2)
plot(Ql,mheatl,'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
xlabel('Heat load [$\frac{J}{cm^2}$]','Interpreter','LaTex','FontSize',14)
ylabel('Forebody heatshield mass[kg]','Interpreter','LaTex','FontSize',14)
grid on