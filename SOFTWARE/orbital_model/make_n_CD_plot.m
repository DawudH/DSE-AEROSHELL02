clc
clear all
close all

constants
addpath('..\matlab2tikz')

fid = fopen('orbit_selection_const_CL.txt','r');
C = textscan(fid,'%f %f %f %f %d %d %d %f');
fclose(fid);

V = C{1,4};
CL = C{1,3};
CD = C{1,2};
max_a = C{1,8};
inorbit = C{1,7};

inorbit_index = find( inorbit == 1 );
V_inorbit = V(inorbit_index);
a_max_inorbit = max_a(inorbit_index);
CD_inorbit = CD(inorbit_index);

different_V = unique(V_inorbit);
different_V = sort(different_V);

figure('name','Different orbits for constant CD')
hold on
grid on
xlabel('$C_D \left[-\right]$','interpreter','LaTeX','fontsize',15)
ylabel('$a_{min} \left[\frac{m}{s^2}\right]$','interpreter','LaTeX','fontsize',18)
%ylim([0 6]*g_earth)
cc = parula(length(different_V) + 1);
legend_str = cell(length(different_V)+1,1);
marker = {'-+'; '-o'; '-*'; '-x'; '-s'; '-d'; '-.'; '-^'; '-v'; '->'; '-<'; '-p'; '-h'};
for i = 1:length(different_V)
        index_V = find( V_inorbit == different_V(i) );
        CD_V = CD_inorbit(index_V);
        [CD_V, index_CD_order] = sort(CD_V);
        a_max_V = a_max_inorbit(index_V);
        a_max_V = a_max_V(index_CD_order);
        
        % remove non min_a terms
        different_CD = unique(CD_V);
        a_plot = zeros(length(different_CD),1);
        
        for k = 1:length(different_CD)
            
            index_CD = find(CD_V == different_CD(k));
            a_max_CD = a_max_V(index_CD);
            a_plot(k) = min(a_max_CD)*g_earth;
            
        end
        
        legend_str{i} = [num2str(different_V(i)) ' m/s'];
        plot(different_CD,a_plot,marker{i},'color',cc(i,:));
end
plot(xlim,[3,3]*g_earth,'-.','color','r');
legend_str{end} = '29.43 m/s^2 (3g)';
legend(legend_str,'location','northoutside');
matlab2tikz('.\LaTeX\n_CD.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
