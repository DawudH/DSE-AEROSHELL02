input = { {'Cd'}, {'qmax'}, {'mf'}, {'M'}, {'T'}, {'Do'}, {'Di'},...
    {'theta'}, {'L'},{'N'},{'I'},{'All material'}};
output ={ {'Decelerator mass [kg]'},'Decelerator mass [-]','Inflation gas mass [-]','Flexible material mass [-]',...
    {'Inflation gas pressure [Pa]','Flexible material mass [kg]','Inflation gas mass [kg]','Ballistic coefficient [$\frac{kg}{m^2}$]'}};
x1 = {'Do'};
x2 = {'qmax'};
y  = {'Ballistic coefficient [$\frac{kg}{m^2}$]'};
plotstyles = { { '2dx1' }, {'2dx2'}, {'3d_all'}, {'3d_stacked'}, {'3d_tension'}, {'3d_trailing'}, {'2d_matlegendx2_all'},{'2d_matlegendx2_stacked'},{'2d_matlegendx2_tension'},...
    {'2d_matlegendx2_trailing'}};
desired_plotstyle = {'3d_mesh'};

%addpath('.\matlab2tikz')
%name = '.\LaTeX\Tikz\stacked_mat.tikz'
figure(1)
Materials = {'Kapton','Kevlar 29','Kevlar 49','M5','Nomex','PBO Zylon','Spectra 2000', 'Technora', 'Upilex-25S','Vectran','Nextel 610'};%,'Aluminium'};
typeoftoroid = [0 1 2];
%x2range = Materials;
%x2range = Materials;
%x2range = typeoftoroid;
%x1range = 0.5:0.1:2.0; %CD
x2range = 1000:500:10000; %dyn pressure
%x1range = 5:2:40; %gas molar mass
%x1range = 1:1:5;
%x1range = 5.01:1:15.01; 
x1range = 5.01:0.5:12;
%x1range = 40:5:80;
%x1range = 2:1:25;
%x2range = 6:1:10;
%x1range = 2:10:300;
%x2range = 273.15:20:573.15;


Cd = 1.5;
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

selected_material_gore = 'Vectran';
selected_material_rad  = 'Vectran';
selected_material_axi  = 'Vectran';
selected_material_gasb = 'Vectran';
selected_material_toro = 'Vectran';

typeoftoroid = 1;

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
            qmax = a; x1lab = {'Peak dynamic pressure [Pa]'};
        elseif strcmp(x1,'mf')
            mf = a; x1lab = {'Regime [-]'};
        elseif strcmp(x1,'M')
            M = a; x1lab = {'Inflation gas molar mass [kg/kmole]'};
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
        elseif strcmp(x1,'Toroid material')
            selected_material_toro = a; x1lab = {'Toroid material [-]'};
        elseif strcmp(x1,'Gas barrier material')
            selected_material_gasb = a; x1lab = {'Gas barrier material [-]'};
        elseif strcmp(x1,'Gore material')
            selected_material_gore = a; x1lab = {'Gore material [-]'};
        elseif strcmp(x1,'Radial strap material')
            selected_material_rad = a; x1lab = {'Radial strap material [-]'};
        elseif strcmp(x1,'Axial strap material')
            selected_material_axi = a; x1lab = {'Axial strap material [-]'};
        elseif strcmp(x1,'Type of toroid')
            typeoftoroid = a; x1lab = {'Type of toroid [-]'};
        elseif strcmp(x1,'All material')
            selected_material_gore = a; selected_material_toro = a ;
            selected_material_gasb = a; selected_material_rad = a  ;
            selected_material_axi = a; x1lab = {'Material [-]'};
        end
    
        if strcmp(x2,'Cd')
            Cd = b; x2lab = {'Drag coefficient [-]'};
        elseif strcmp(x2,'qmax')
            qmax = b; x2lab = {'Peak dynamic pressure [Pa]'};
        elseif strcmp(x2,'mf')
            mf = b; x2lab = {'Regime [-]'};
        elseif strcmp(x2,'M')
            M = b; x2lab = {'Inflation gas molar mass [kg/kmole]'};
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
        elseif strcmp(x2,'Toroid material')
            selected_material_toro = b; x2lab = {'Toroid material [-]'};
        elseif strcmp(x2,'Gas barrier material')
            selected_material_gasb = b; x2lab = {'Gas barrier material [-]'};         
        elseif strcmp(x2,'Gore material')
            selected_material_gore = b; x2lab = {'Gore material [-]'};
        elseif strcmp(x2,'Radial strap material')
            selected_material_rad = b; x2lab = {'Radial strap material [-]'};
        elseif strcmp(x2,'Axial strap material')
            selected_material_axi = b; x2lab = {'Axial strap material [-]'};
        elseif strcmp(x2,'Type of toroid')
            typeoftoroid = b; x2lab = {'Type of toroid [-]'};
        elseif strcmp(x2,'All material')
            selected_material_gore = b; selected_material_toro = b ;
            selected_material_gasb = b; selected_material_rad = b  ;
            selected_material_axi = b; x2lab = {'Material [-]'};
        
        end
        
        [mtoroidbar,mgasbar,mfactor,mtotalbar,mtotal,pmin,mfilmbar,mgasbarrierbar,mcoatingbar,mgoresbar,maxialbar,mradialbar,minflsysbar] = ST_TC_TIAD_massfunction(Cd,qmax,mf,M,T,...
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


% if strcmp(x1,{'Cd'}) && not(strcmp(x2,{'Do'}))     
%     CdA = x1range.*pi/4*Do^2;
%     CdA = repmat(CdA',[1 length(x2range)]);
%     BC_1 = m_total_1./(CdA);BC_2 = m_total_2./(CdA);BC_3 = m_total_3./(CdA);
% elseif strcmp(x2,{'Cd'}) && not(strcmp(x1,{'Do'})) 
%     CdA = x2range.*pi/4*Do^2;
%     CdA = repmat(CdA',[1 length(x1range)]);
%     BC_1 = m_total_1./(CdA);BC_2 = m_total_2./(CdA);BC_3 = m_total_3./(CdA);
% end

if strcmp(y,'Decelerator mass [kg]')
    y1 = m_total_1; y2 = m_total_2 ; y3 = m_total_3;
elseif strcmp(y,'Decelerator mass [-]')
    y1 = m_totalbar_1; y2 = m_totalbar_2 ; y3 = m_totalbar_3;
elseif strcmp(y,'Decelerator mass [-]')
    y1 = m_totalbar_1; y2 = m_totalbar_2 ; y3 = m_totalbar_3;   
elseif strcmp(y,'Inflation gas mass [-]')
    y1 = m_gasbar_1  ; y2 = m_gasbar_2   ; y3 = m_gasbar_3  ;
elseif strcmp(y,'Inflation gas mass [kg]')
    y1 = m_gasbar_1.*m_total_1./m_totalbar_1 ; y2 = m_gasbar_2.*m_total_2./m_totalbar_2  ; y3 = m_gasbar_3.*m_total_3./m_totalbar_3 ;
elseif strcmp(y,'Flexible material mass [-]')
    y1 = m_flexbar_1 ; y2 = m_flexbar_2  ; y3 = m_flexbar_3 ;
elseif strcmp(y,'Flexible material mass [kg]')
    y1 = m_flexbar_1.*m_total_1./m_totalbar_1 ; y2 = m_flexbar_2.*m_total_2./m_totalbar_2  ; y3 = m_flexbar_3.*m_total_3./m_totalbar_3 ;
elseif strcmp(y,'Inflation gas pressure [Pa]')
    y1 = pmin_1      ; y2 = pmin_2       ; y3 = pmin_3      ;
elseif strcmp(y,'Ballistic coefficient [$\frac{kg}{m^2}$]')
    CdA = Cd.*pi/4.*x1range.^2;
    CdA = repmat(CdA',[1 length(x2range)]);
    y1 = m_total_1./CdA; y2 = m_total_2./CdA; y3 = m_total_3./CdA;
end

if strcmp(desired_plotstyle,'3d_all')   
    rot_angle = -50;
    subplot(131)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material'))  || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y1');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y1');
    end
    colormap parula
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    subplot(132)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material'))   || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y2');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y2');
    end
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    subplot(133)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y3');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y3');
    end
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    
    h = colorbar;
    set(h,'Position',[0.07 0.1 0.02 0.8]);    
    set(h, 'YAxisLocation','left')
    ylabel(h,y,'Interpreter','LaTex','FontSize',13);
  
elseif strcmp(desired_plotstyle,'3d_stacked')   
    rot_angle = -50;
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y1');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y1');
    end
    colormap parula
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    h = colorbar;
    %set(h,'Position',[0.07 0.1 0.02 0.8]);    
    %set(h, 'YAxisLocation','left')
    ylabel(h,y,'Interpreter','LaTex','FontSize',13);

elseif strcmp(desired_plotstyle,'3d_tension')   
    rot_angle = -50;
    
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material')) 
        C = contourf(1:length(x1range),1:length(x2range),y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y2');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y2');
    end
    
    colormap parula
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    h = colorbar;
    %set(h,'Position',[0.07 0.1 0.02 0.8]);    
    %set(h, 'YAxisLocation','left')
    ylabel(h,y,'Interpreter','LaTex','FontSize',13);

elseif strcmp(desired_plotstyle,'3d_trailing')   
    rot_angle = -50;
    
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y3');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = contourf(x1range,x2range,y3');
    end
    colormap parula
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    h = colorbar;
    %set(h,'Position',[0.07 0.1 0.02 0.8]);    
    %set(h, 'YAxisLocation','left')
    ylabel(h,y,'Interpreter','LaTex','FontSize',13);

      
elseif strcmp(desired_plotstyle,'2dx1')
    if strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        plot(1:length(x1range),y1(1:end,1)','d','Color','k','LineStyle','-','LineWidth',1);
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        hold on
        if strcmp(x1,{'N'})==0
            plot(1:length(x1range),y2(1:end,1)','s','Color','m','LineStyle','-','LineWidth',1)
            hold on
            plot(1:length(x1range),y3(1:end,1)','v','Color','b','LineStyle','-','LineWidth',1)
        end
    else
        plot(x1range,y1(1:end,1)','d','Color','k','LineStyle','-','LineWidth',1)
        hold on
         if strcmp(x1,{'N'})==0
            plot(x1range,y2(1:end,1)','s','Color','m','LineStyle','-','LineWidth',1)
            hold on
            plot(x1range,y3(1:end,1)','v','Color','b','LineStyle','-','LineWidth',1)
            hold on
            %plot(Dl,BCl(end:-1:1),'*','Color','r','LineStyle','-','LineWidth',1)
        end
    end
    legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
elseif strcmp(desired_plotstyle,'2dx2')
    if strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        plot(1:length(x2range),y1(1,1:end)','d','Color','k','LineStyle','-','LineWidth',1);
        set(gca, 'XTick',1:length(x2range), 'XTickLabel',x2range) 
        hold on
        if strcmp(x2,{'N'})==0
            plot(1:length(x2range),y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
            hold on
            plot(1:length(x2range),y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
        end
    else
        plot(x2range,y1(1,1:end)','d','Color','k','LineStyle','-','LineWidth',1)
        hold on
        if strcmp(x2,{'N'})==0
            plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
            hold on
            plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
        end
    end
    legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    xlabel(x2lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
    
    %WIP
elseif strcmp(desired_plotstyle,'2d_matlegendx2_all')
    
    subplot(131)
    markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
    %colors = {'k','b','r','g','y','c','m'[.5 .6 .7],[.8 .2 .6]};
    for mat = 1:length(x2range)
        plot(x1range,y1(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    end
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
    
    subplot(132)
    for mat = 1:length(x2range)
        plot(x1range,y2(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    %if strcmp(x2,{'N'})==0
    %        plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
    %        hold on
    %        plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
    %legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
    end
    
    subplot(133)
    for mat = 1:length(x2range)
        plot(x1range,y3(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    %if strcmp(x2,{'N'})==0
    %        plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
    %        hold on
    %        plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
    %legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    end
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    legend(Materials,'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
    
    
    
elseif strcmp(desired_plotstyle,'2d_matlegendx2_stacked')
    
    markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
    %colors = {'k','b','r','g','y','c','m'[.5 .6 .7],[.8 .2 .6]};
    for mat = 1:length(x2range)
        plot(x1range,y1(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    end
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    legend(Materials,'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on
    
elseif strcmp(desired_plotstyle,'2d_matlegendx2_tension')
    
    markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
    %colors = {'k','b','r','g','y','c','m'[.5 .6 .7],[.8 .2 .6]};
  
    for mat = 1:length(x2range)
        plot(x1range,y2(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    end
    %if strcmp(x2,{'N'})==0
    %        plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
    %        hold on
    %        plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
    %legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    legend(Materials,'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on

elseif strcmp(desired_plotstyle,'2d_matlegendx2_trailing')
    
    markers = {'+','o','*','h','x','s','d','^','v','>','<','p','.'};
    %colors = {'k','b','r','g','y','c','m'[.5 .6 .7],[.8 .2 .6]};
    
    for mat = 1:length(x2range)
        plot(x1range,y3(1:end,mat)',markers{mod(mat,numel(markers))+1},'Color','k','Linestyle','-','LineWidth',1.1,'MarkerSize',6.5)
        hold on
    %if strcmp(x2,{'N'})==0
    %        plot(x2range,y2(1,1:end)','s','Color','m','LineStyle','-','LineWidth',1)
    %        hold on
    %        plot(x2range,y3(1,1:end)','v','Color','b','LineStyle','-','LineWidth',1)
    %legend({'Stacked toroid','Tension cone', 'Trailing IAD'},'Interpreter','LaTex','FontSize',12,2)
    end
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    legend(Materials,'Interpreter','LaTex','FontSize',12,2)
    xlabel(x1lab,'Interpreter','LaTex','FontSize',14)
    ylabel(y,'Interpreter','LaTex','FontSize',14)
    grid on

elseif strcmp(desired_plotstyle,'3d_mesh')   
    rot_angle = -50;
    %subplot(131)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y1');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material'))  || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y1');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = surface(x1range,x2range,y1');
        hold on
    end
    colormap hsv
    title({'\textbf{Stacked toroid}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    %subplot(132)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material'))   || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y2');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y2');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = surface(x1range,x2range,y2');
        hold on
    end
    colormap parula
    title({'\textbf{Tension cone}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    %subplot(133)
    if (strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material'))) && strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))  || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),1:length(x2range),y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    elseif strcmp(x1,'Toroid material') || (strcmp(x1,'Gas barrier material')) || (strcmp(x1,'Axial strap material')) || (strcmp(x1,'Radial strap material')) || (strcmp(x1,'Gore material')) || (strcmp(x1,'All material'))
        C = contourf(1:length(x1range),x2range,y3');
        set(gca, 'XTick',1:length(x1range), 'XTickLabel',x1range) 
    elseif strcmp(x2,'Toroid material') || (strcmp(x2,'Gas barrier material')) || (strcmp(x2,'Axial strap material')) || (strcmp(x2,'Radial strap material')) || (strcmp(x2,'Gore material')) || (strcmp(x2,'All material'))
        C = contourf(x1range,1:length(x2range),y3');
        set(gca, 'YTick',1:length(x2range), 'YTickLabel',x2range) 
    else
        C = surface(x1range,x2range,y3');
    end
    title({'\textbf{Trailing IAD}'},'Interpreter','LaTex','FontSize',13)
    hxlabel=xlabel(x1lab,'Interpreter','LaTex','FontSize',13);
    hylabel=ylabel(x2lab,'Interpreter','LaTex','FontSize',13);
    hzlabel=zlabel(y,'Interpreter','LaTex','FontSize',13);
    %set(hylabel,'rotation',rot_angle)
    %set(hxlabel,'rotation',90+rot_angle)
    grid on
    
    legend('Stacked','Tension cone','Trailing IAD')
    
    
    %h = colorbar;
    %set(h,'Position',[0.07 0.1 0.02 0.8]);    
    %set(h, 'YAxisLocation','left')
    %ylabel(h,y,'Interpreter','LaTex','FontSize',13);

end
%matlab2tikz(name,'height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
    