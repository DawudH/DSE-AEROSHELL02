clear
clc
%close all

constants

load('results.mat');
figure('name','under 3g')
xlabel('Flight path angle')
ylabel('Angle of attack')
hold on
for k = 1:length(results.A)
    if results.A(k)<3*g_earth
        plot(results.GAMMA(k),results.ALPHA(k)*180/pi,'bx')
    else
         plot(results.GAMMA(k),results.ALPHA(k)*180/pi,'rx')
    end
end