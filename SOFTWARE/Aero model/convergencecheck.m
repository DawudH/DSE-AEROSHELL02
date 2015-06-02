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
shapetexts.ballute = 'ballute';

a = 300;
gamma = 1.4;
center = zeros(3,1);
rho = 1e-3;
T = 150;
qarray = [5:2:101];
qmaxarray = zeros(size(qarray));
cd = zeros(size(qarray));
cl = zeros(size(qarray));

for i = 1:length(qarray)
    qarray(i)
    
    [ TriGeom, A, center ] = generategeometry( qarray(i) );
    geom = aeroGeometry(TriGeom, A);
    
    mod = modnewtonian(geom, gamma, a, center, rho, T);
    mod = mod.calcAeroangle(7000,deg2rad(20),deg2rad(0));

    Tw = 500*ones(size(mod.Cpdist_array(:,end)));
    [Tmax, qmax, qw] = mod.calcStagnationHeatFlux(Tw);
    qmaxarray(i) = qmax;
    cd(i) = mod.CRA_aero_array(1,end);
    cl(i) = mod.CRA_aero_array(3,end);
    
end

figure;
plot(qarray, qmaxarray);
ylabel('qmax');

figure;
plot(qarray, cl);
ylabel('cl');
figure;
plot(qarray, cd);
ylabel('cd');


disp('q:')
ratio = qmaxarray/qmaxarray(end);
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