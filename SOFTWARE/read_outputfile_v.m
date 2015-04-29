% read file and plot
clc
clear all
close all

fid = fopen('orbit_true_or_false_v.txt','r');
C = textscan(fid,'%11.1f %4.2f %d %d %d %f %f');
fclose(fid);


figure('name','results')
hold on
grid on
inorbit = C{5};
rx = C{1};
CD = C{2};
delta_v = C{7};
acc = C{6}

for i=1:length(delta_v)
   %delta_v(i)
    if acc(i)<3
        scatter3(rx(i),CD(i),delta_v(i),'*')
    end
    
end