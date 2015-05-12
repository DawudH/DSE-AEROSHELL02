% close all;
horizontalplate = 'horizontalplate';
verticalplate = 'verticalplate';
sphere12m = 'sphere12m';
pastille12m3m = 'pastille12m1.5m';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
q = 10;

alpha0 = 0; %degrees
dalpha = 1; %degrees
alphaend = 15; %degrees
qarray = [4 5 6 7 8 9 10,10:1:40, 42:2:60, 65:5:100];
CD = [];

for q = qarray
    q
    [ coords, tri, A ] = generategeometry( 'deg30cone', q );

    
    mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);
    % mod = mod.alphasweep(7e3, 0, deg2rad(alpha0), deg2rad(alphaend), deg2rad(dalpha));

    mod = mod.calcAeroangle(7e3,deg2rad(0),0);
    CD = [CD, mod.CR_aero_array(1,1)];
end
CD;
    % mod.plotCp(true, true);
    % mod.CR_aero_array
figure;
plot(qarray, CD);
% mod.plots(rad2deg(mod.alpha_array), 'alpha', {{'cl'}, {'cd'}, {'clcd'}, {'cx'}, {'cz'},{'cmy'},{'q'}, {'T'}});

ratio = CD/CD(end);
disp(strcat('Best q for this geometry: ', qarray(find(ratio(ratio>0.999))));
% for i = 2:length(CD)
% i;
% CD(i)/CD(end);
% end
    