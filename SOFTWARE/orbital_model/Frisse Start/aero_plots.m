clear
close all
clc

% generate cl alpha plots
added_paths
constants

alpha = -25:1:25;
alpha = alpha * pi / 180;

torus = aeroProperties('torus');
isotensoid = aeroProperties('isotensoid');
irve = aeroProperties('irve');
apollo = aeroProperties('apollo');
pastille = aeroProperties('pastille');

[o_torus.CLA, o_torus.CDA, o_torus.CMYA] = torus.aeroCoeffs(alpha);
[o_isotensoid.CLA, o_isotensoid.CDA, o_isotensoid.CMYA] = isotensoid.aeroCoeffs(alpha);
[o_irve.CLA, o_irve.CDA, o_irve.CMYA] = irve.aeroCoeffs(alpha);
[o_apollo.CLA, o_apollo.CDA, o_apollo.CMYA] = apollo.aeroCoeffs(alpha);
[o_pastille.CLA, o_pastille.CDA, o_pastille.CMYA] = pastille.aeroCoeffs(alpha);

marker = {'-d'; '-o'; '-*'; '-x'; '-s'; '-+'; '-.'; '-^'; '-v'; '->'; '-<'; '-p'; '-h'};
cc = parula(8);
figure('name','aero coeffs')
hold on
grid on
alpha = alpha * 180/pi;
plot(alpha,o_torus.CLA / S,marker{1},'color',cc(1,:)) 
plot(alpha,o_isotensoid.CLA / S,marker{2},'color',cc(2,:)) 
plot(alpha,o_irve.CLA / S,marker{3},'color',cc(3,:)) 
plot(alpha,o_apollo.CLA / S,marker{4},'color',cc(4,:)) 
plot(alpha,o_pastille.CLA / S,marker{5},'color',cc(5,:)) 
legend('torus','isotensoid','irve','apollo','pastille')
