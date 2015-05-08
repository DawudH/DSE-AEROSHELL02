% read file and plot
clc
clear all
close all

constants
addpath('..\matlab2tikz')

fid = fopen('orbit_selection_range.txt','r');
C = textscan(fid,'%f %f %f %f %d %d %d %f');
fclose(fid);

RX = C{1,1};
V = C{1,4};
CL = C{1,3};
CD = C{1,2};
max_a = C{1,8};
inorbit = C{1,7};

inorbit_index = find( inorbit == 1 );
V_inorbit = V(inorbit_index);
a_max_inorbit = max_a(inorbit_index);
CL_inorbit = CL(inorbit_index);
RX_inorbit = RX(inorbit_index);

distance = max(RX_inorbit) - min(RX_inorbit); 


fig = figure('name','results');
hold on
grid on
v_text = 3.2*g_earth;
text(min(RX_inorbit),v_text,'$\leftarrow$','HorizontalAlignment','left','Interpreter','LaTeX','fontsize',15)
text(max(RX_inorbit),v_text,'$\rightarrow$','HorizontalAlignment','right','Interpreter','LaTeX','fontsize',15)
text(min(RX_inorbit) + distance/2,v_text,[' $' num2str(distance) '\, \left[m\right]$'],'HorizontalAlignment','center','Interpreter','LaTeX','fontsize',15)
ylim([0 6]*g_earth)
xlim([min(RX_inorbit)-distance/2 max(RX_inorbit)+distance/2])
cc = parula(5);
plot(RX_inorbit,a_max_inorbit*g_earth,'-o','color',cc(2,:))
xlabel('$R_x \left[-\right]$','interpreter','LaTeX','fontsize',15)
ylabel('$a \left[\frac{m}{s^2}\right]$','interpreter','LaTeX','fontsize',18)
matlab2tikz('.\LaTeX\range_plot.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);


