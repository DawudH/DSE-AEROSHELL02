% read file and plot
clc
clear all
close all

fid = fopen('orbit_true_or_false_L.txt','r');
C = textscan(fid,'%11.1f %4.2f %d %d %d %f');
fclose(fid);


fig = figure('name','results')
hold on
grid on
inorbit = C{5};
rx = C{1};
CD = C{2};
a = C{6};

k = 1;
for i=1:length(inorbit)
   inorbit(i)
    if (inorbit(i) == 1)
%         plot(rx(i),CD(i),'*')
        rx_orb(k) = rx(i);
        CD_orb(k) = CD(i);
        a_orb(k) = a(i);
        k = k+1;
    end
    
end
[Xm,Ym] = meshgrid(rx_orb,CD_orb);
[X, Y,Z] = griddata(rx_orb,CD_orb,a_orb,Xm,Ym);
plot3(X,Y,Z,'o')
%colormap(fig,jet)
