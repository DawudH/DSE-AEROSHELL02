% read file and plot
clc
clear all
close all

fid = fopen('orbit_true_or_false.txt','r');
C = textscan(fid,'%11.1f %4.2f %d %d %d %f');
fclose(fid);


figure('name','results')
hold on
grid on
inorbit = C{5};
rx = C{1};
CD = C{2};

for i=1:length(inorbit)
   inorbit(i)
    if (inorbit(i) == 1)
        plot(rx(i),CD(i),'*')
    end
    
end