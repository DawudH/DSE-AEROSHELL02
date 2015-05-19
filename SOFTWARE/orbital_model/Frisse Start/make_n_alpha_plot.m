clc
clear all
close all

constants
addpath('..\..\matlab2tikz')

fid = fopen('orbit_alpha_fine.txt','r');
C = textscan(fid,'%f %f %f %f %f %d %d %d %f');
fclose(fid);
% filestr{1} = sprintf('rx(m) \t ry(m) \t v \t dt  \t alpha(deg) \t  crash \t flyby \t orbit \t max_accel \n');
V = C{1,3};
alpha = C{1,5};
max_a = C{1,9};
inorbit = C{1,8};

inorbit_index = find( inorbit == 1 );
V_inorbit = V(inorbit_index);
a_max_inorbit = max_a(inorbit_index);
alpha_inorbit = alpha(inorbit_index);

different_V = unique(V_inorbit);

figure('name','Different orbits for certain alpha without control')
hold on
grid on
xlabel('$\alpha \left[^\circ\right]$','interpreter','LaTeX','fontsize',15)
ylabel('$a_{max} \left[\frac{m}{s^2}\right]$','interpreter','LaTeX','fontsize',15)
%ylim([0 6]*g_earth)
cc = parula(length(different_V) + 3);
legend_str = cell(length(different_V)*3+1,1);
marker = {'-d'; '-o'; '-*'; '-x'; '-s'; '-+'; '-.'; '-^'; '-v'; '->'; '-<'; '-p'; '-h'};
j = 1;
for i = 1:length(different_V)
        index_V = find( V_inorbit == different_V(i) );
        alpha_V = alpha_inorbit(index_V);
        [alpha_V, index_alpha_order] = sort(alpha_V);
        a_max_V = a_max_inorbit(index_V);
        a_max_V = a_max_V(index_alpha_order);
        
        % remove non min_a terms
        different_alpha = unique(alpha_V);
        a_plot_min = zeros(length(different_alpha),1);
        a_plot_max = zeros(length(different_alpha),1);
        
        for k = 1:length(different_alpha)
            
            index_alpha = find(alpha_V == different_alpha(k));
            a_max_alpha = a_max_V(index_alpha);
            a_plot_min(k) = min(a_max_alpha)*g_earth;
            a_plot_max(k) = max(a_max_alpha)*g_earth;
            
        end
        
        legend_str{j} = [num2str(different_V(i)) ' [m/s], lower'];
        legend_str{j+1} = [num2str(different_V(i)) ' [m/s], upper'];
        legend_str{j+2} = [num2str(different_V(i)) ' [m/s], in orbit area'];
        j = j + 3;
        plot(different_alpha,a_plot_min,marker{i},'color',cc(i,:));
        plot(different_alpha,a_plot_max,marker{i+4},'color',cc(i,:));
        fill([different_alpha; flip(different_alpha)],[a_plot_min; flip(a_plot_max)],cc(i,:),'FaceAlpha',0.2)
        %fill([different_alpha; flip(different_alpha)],[a_plot_min; flip(a_plot_max)],cc(i,:))
        
end
plot(xlim,[3,3]*g_earth,'-.','color',cc(i+1,:));
legend_str{end} = '29.43 m/s^2 (3g)';
legend(legend_str,'location','northeast');
matlab2tikz('.\LaTeX\n_alpha.tikz','height','\figureheight','width','\figurewidth','showInfo', false,'checkForUpdates',false);
