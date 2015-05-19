% close all;

shapetexts.horizontalplate = 'horizontalplate';
shapetexts.verticalplate = 'verticalplate';
shapetexts.sphere12m = 'sphere12m';
shapetexts.pastille12m15m = 'pastille12m1.5m';
shapetexts.deg60cone = 'deg60cone';
shapetexts.deg30cone = 'deg30cone';
shapetexts.irvevalidation = 'irvevalidation';
shapetexts.apollovalidation = 'apollovalidation';
shapetexts.torus = 'torus';
shapetexts.concept_irve = 'concept_irve';
shapetexts.concept_apollo = 'concept_apollo';
shapetexts.concept_isotensoid = 'concept_isotensoid';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
qarray = [5:2:43,47:2:81];
qmax = zeros(size(qarray));
cd = zeros(size(qarray));
cl = zeros(size(qarray));

for i = 1:length(qarray)
    q = qarray(i);
    disp(strcat('Current q:', num2str(q)));
    [ coords, tri, A ] = generategeometry( shapetexts.concept_isotensoid, q );
    
    mod = modnewtonian( coords, tri, gamma, a, center, rho, T, A);

    mod = mod.calcAeroangle(7e3,deg2rad(1),0);
    qmax(i) = mod.qmax_array(1,1);
    cd(i) = mod.CRA_aero_array(1,1);
    cl(i) = mod.CRA_aero_array(3,1);
end

figure;
plot(qarray, qmax);

disp('q:')
ratio = qmax/qmax(end);
bestarray = find((1.001>ratio) & (ratio>0.999));
disp(strcat('q required for 99.9% convergence: ', num2str(qarray(bestarray(1)))));
bestarray = find((1.005>ratio) & (ratio>0.995));
disp(strcat('q required for 99.5% convergence: ', num2str(qarray(bestarray(1)))));
bestarray = find((1.01>ratio) & (ratio>0.99));
disp(strcat('q required for 99% convergence: ', num2str(qarray(bestarray(1)))));

disp('');
disp('cl:');
ratio = cl/cl(end);
bestarray = find((1.001>ratio) & (ratio>0.999));
disp(strcat('q required for 99.9% convergence: ', num2str(qarray(bestarray(1)))));
bestarray = find((1.005>ratio) & (ratio>0.995));
disp(strcat('q required for 99.5% convergence: ', num2str(qarray(bestarray(1)))));  
bestarray = find((1.01>ratio) & (ratio>0.99));
disp(strcat('q required for 99% convergence: ', num2str(qarray(bestarray(1)))));

disp('');
disp('cd:');
ratio = cd/cd(end);
bestarray = find((1.001>ratio) & (ratio>0.999));
disp(strcat('q required for 99.9% convergence: ', num2str(qarray(bestarray(1)))));
bestarray = find((1.005>ratio) & (ratio>0.995));
disp(strcat('q required for 99.5% convergence: ', num2str(qarray(bestarray(1)))));  
bestarray = find((1.01>ratio) & (ratio>0.99));
disp(strcat('q required for 99% convergence: ', num2str(qarray(bestarray(1)))));