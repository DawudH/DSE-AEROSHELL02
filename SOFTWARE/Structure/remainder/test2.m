input = { {'Cd'}, {'qmax'}, {'mf'}, {'M'}, {'T'}, {'Do'}, {'Di'},...
    {'theta'}, {'L'},{'N'},{'I'}};
output ={ {'Total mass [kg]'},'Total mass [-]','Inflation gas mass [-]','Flexible material mass [-]',...
    {'Inflation gas pressure [Pa]'}};
x1 = {'qmax'};
x2 = {'theta'};
y  = {'Total mass [kg]'};
desired_plotstyle = {'2dx2'};

x1range = 0.5:0.1:2.0;
x1range = 1000:100:3000;
%x2range = 5.01:0.5:12;
x2range = 40:5:80;


Cd = 0.5;
qmax = 3000;
mf = 2;
M = 22;
T = 273.15;
Do= 12;
Di= 5;
theta = 60;
L = 1;
N = 8;
I = 18;

selected_material_gore = 'Kevlar 49';
selected_material_rad  = 'Kevlar 49';
selected_material_axi  = 'Kevlar 49';
selected_material_gasb = 'Kevlar 49';
selected_material_toro = 'Kevlar 49';

typeoftoroid = 0;

%x1range = Cdrange ; x2range = Dorange ;


m_total_1 = zeros(length(x1range),length(x2range));
m_total_2 = zeros(length(x1range),length(x2range));
m_total_3 = zeros(length(x1range),length(x2range));

m_totalbar_1 = zeros(length(x1range),length(x2range));
m_totalbar_2 = zeros(length(x1range),length(x2range));
m_totalbar_3 = zeros(length(x1range),length(x2range));

m_gasbar_1 = zeros(length(x1range),length(x2range));
m_gasbar_2 = zeros(length(x1range),length(x2range));
m_gasbar_3 = zeros(length(x1range),length(x2range));

m_flexbar_1 = zeros(length(x1range),length(x2range));
m_flexbar_2 = zeros(length(x1range),length(x2range));
m_flexbar_3 = zeros(length(x1range),length(x2range));

pmin_1 = zeros(length(x1range),length(x2range));
pmin_2 = zeros(length(x1range),length(x2range));
pmin_3 = zeros(length(x1range),length(x2range));

i = 1;
while i<=length(x1range)
    j = 1;
    a = x1range(i);
    while j<=length(x2range)
        b = x2range(j);
        
        if strcmp(x1,'Cd')
            Cd = a; x1lab = {'Drag coefficient [-]'};
        elseif strcmp(x1,'qmax')
            qmax = a; x1lab = {'Dynamic pressure [Pa]'};
        elseif strcmp(x1,'mf')
            mf = a; x1lab = {'Regime [-]'};
        elseif strcmp(x1,'M')
            M = a; x1lab = {'Inflation gas molar mass [\frac{kg}{kmole}]'};
        elseif strcmp(x1,'T')
            T = a; x1lab = {'Inflation gas temperature [K]'};
        elseif strcmp(x1,'Do')
            Do= a; x1lab = {'Deployed diameter [m]'};
        elseif strcmp(x1,'Di')
            Di= a; x1lab = {'Centerbody diameter [m]'};
        elseif strcmp(x1,'theta')
            theta=a; x1lab = {'Half-cone angle [deg]'};
        elseif strcmp(x1,'L')
            L = a; x1lab = {'Characteristic length [m]'};
        elseif strcmp(x1,'N')
            N = a; x1lab = {'Number of toroids [-]'};
        elseif strcmp(x1,'I')
            I = a; x1lab = {'Number of radial straps [-]'};
        end
        
        if strcmp(x2,'Cd')
            Cd = b; x2lab = {'Drag coefficient [-]'};
        elseif strcmp(x2,'qmax')
            qmax = b; x2lab = {'Dynamic pressure [Pa]'};
        elseif strcmp(x2,'mf')
            mf = b; x2lab = {'Regime [-]'};
        elseif strcmp(x2,'M')
            M = b; x2lab = {'Inflation gas molar mass [\frac{kg}{kmole}]'};
        elseif strcmp(x2,'T')
            T = b; x2lab = {'Inflation gas temperature [K]'};
        elseif strcmp(x2,'Do')
            Do= b; x2lab = {'Deployed diameter [m]'};
        elseif strcmp(x2,'Di')
            Di= b; x2lab = {'Centerbody diameter [m]'};
        elseif strcmp(x2,'theta')
            theta=b; x2lab = {'Half-cone angle [deg]'};
        elseif strcmp(x2,'L')
            L = b; x2lab = {'Characteristic length [m]'};
        elseif strcmp(x2,'N')
            N = b; x2lab = {'Number of toroids [-]'};
        elseif strcmp(x2,'I')
            I = b; x2lab = {'Number of radial straps [-]'};
        end
        
        [mtoroidbar,mgasbar,mfactor,mtotalbar,mtotal,pmin] = ST_TC_TIAD_massfunction(Cd,qmax,mf,M,T,...
Do,Di,theta,L,N,I,...
selected_material_gore,selected_material_rad,selected_material_axi,...
selected_material_gasb,selected_material_toro,typeoftoroid);
        
        m_total_1(i,j)=mtotal(1);
        m_total_2(i,j)=mtotal(2);
        m_total_3(i,j)=mtotal(3);
        
        m_totalbar_1(i,j)=mtotalbar(1);
        m_totalbar_2(i,j)=mtotalbar(2);
        m_totalbar_3(i,j)=mtotalbar(3);
        
        m_gasbar_1(i,j)=mgasbar(1);
        m_gasbar_2(i,j)=mgasbar(2);
        m_gasbar_3(i,j)=mgasbar(3);
        
        m_flexbar_1(i,j)=mtoroidbar(1);
        m_flexbar_2(i,j)=mtoroidbar(2);
        m_flexbar_3(i,j)=mtoroidbar(3);
        
        pmin_1(i,j)=pmin(1);
        pmin_2(i,j)=pmin(2);
        pmin_3(i,j)=pmin(3);
        
        j = j + 1;
    end
    
    i = i + 1;
end

% if strcmp(x2,'Cd')
%             Cd = b;
%        elseif strcmp(x2,'qmax')
%             qmax = b;
%         elseif strcmp(x2,'mf')
%             mf = b;
%         elseif strcmp(x2,'M')
%             M = b;
%         elseif strcmp(x2,'T')
%             T = b;
%         elseif strcmp(x2,'Do')
%             Do= b;
%         elseif strcmp(x2,'Di')
%             Di= b;
%         elseif strcmp(x2,'theta')
%             theta=b;
%         elseif strcmp(x2,'L')
%             L = b;
%         elseif strcmp(x2,'N')
%             N = b;
%         elseif strcmp(x2,'I')
%             I = b;
% end

if strcmp(y,'Total mass [kg]')
    y1 = m_total_1; y2 = m_total_2 ; y3 = m_total_3;
elseif strcmp(y,'Total mass [-]')
    y1 = m_totalbar_1; y2 = m_totalbar_2 ; y3 = m_totalbar_3;
elseif strcmp(y,'Total mass [-]')
    y1 = m_totalbar_1; y2 = m_totalbar_2 ; y3 = m_totalbar_3;   
elseif strcmp(y,'Inflation gas mass [-]')
    y1 = m_gasbar_1  ; y2 = m_gasbar_2   ; y3 = m_gasbar_3  ;
elseif strcmp(y,'Flexible material mass [-]')
    y1 = m_flexbar_1 ; y2 = m_flexbar_2  ; y3 = m_flexbar_3 ;
elseif strcmp(y,'Inflation gas pressure [Pa]')
    y1 = pmin_1      ; y2 = pmin_2       ; y3 = pmin_3      ;
end

if strcmp(desired_plotstyle,'3d')   
    rot_angle = -50;
    subplot(131)
    surf(x1range,x2range,y1')
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex');
    hylabel=ylabel(x2lab,'Interpreter','LaTex');
    hzlabel=zlabel(y,'Interpreter','LaTex');
    set(hylabel,'rotation',rot_angle)
    set(hxlabel,'rotation',90+rot_angle)
    grid on
    subplot(132)
    surf(x1range,x2range,y2')
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex');
    hylabel=ylabel(x2lab,'Interpreter','LaTex');
    hzlabel=zlabel(y,'Interpreter','LaTex');
    set(hylabel,'rotation',rot_angle)
    set(hxlabel,'rotation',90+rot_angle)
    grid on
    subplot(133)
    surf(x1range,x2range,y3')
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex');
    hylabel=ylabel(x2lab,'Interpreter','LaTex');
    hzlabel=zlabel(y,'Interpreter','LaTex');
    set(hylabel,'rotation',rot_angle)
    set(hxlabel,'rotation',90+rot_angle)
    grid on
  
elseif strcmp(desired_plotstyle,'2dx1')
    plot(x1range,y1(1:end,1)','d','Color','k','LineStyle','-','LineWidth',1)
    hold on
    plot(x1range,y2(1:end,1)','s','Color','m','LineStyle','-','LineWidth',1)
    hold on
    plot(x1range,y3(1:end,1)','v','Color','b','LineStyle','-','LineWidth',1)
    legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
elseif strcmp(desired_plotstyle,'2dx2')
    plot(x2range,y1(1,1:end)','d','Color','k','LineStyle','-','LineWidth',1)
    hold on
    plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
    hold on
    plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
    legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    xlabel(x2lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
end
   
    