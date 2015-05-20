clear all
close all
clc

%% Determination of heat loads

Tatm = cell(5,1);
qearo = cell(5,1);


% load the concept atmosphere temps and the aerodynamic heat fluxes 
load('heatflux_ballute.mat');
qearo{1,1} = qmax;
Tatm{1,1} = T;
load('heatflux_isotensoid.mat');
qearo{2,1} = qmax;
Tatm{2,1}   = T;
load('heatflux_rigid.mat');
qearo{3,1} = qmax;
Tatm{3,1}    = T;
load('heatflux_stackedtoroid_tensioncone.mat');
qearo{4,1} = qmax;
Tatm{4,1}   = T;
qearo{5,1} = qmax;
Tatm{5,1}    = T;

%Determine dt
dt = t(length(t))/length(t);

% plot the temp
T = Tatm{1,1};
% figure(1)
% plot(t,T)
% hold off



% Define different wall temperatures
min  = 500;
mout = 1500;

for m = min:mout;

   dT{m-(min-1),1} = m.^4 - T.^4;
   
end
 

% dT      = cell(5,1);
% dT{1,1} = T.^4 - T.^4;
% dT{2,1} = 0500^4 - T.^4;
% dT{3,1} = 1000^4 - T.^4;
% dT{4,1} = 1500^4 - T.^4;
% dT{5,1} = 2000^4 - T.^4;

% Material emissivity per concept
eps      = zeros(5,1);
eps(1,1) = 0.443;
eps(2,1) = 0.443;
eps(3,1) = 0.443;
eps(4,1) = 0.443;
eps(5,1) = 0.443;

% Radiation determination
sig = 5.670373e-8;

% Define var
qoutc = cell(5,1);
qtot = cell(5,1);
Qtot = cell(5,1);
qtoto = zeros(length(t),1);
Qmax = zeros(5,1);

% sensativity for different concepts
for k = 1:length(dT)
    
    % for loop, to analyse each concept
    for i = 1:length(Tatm)  
       %Detrmine the output flux 
       qoutc{i,1} = (1/10000).*eps(i,1)*sig.*dT{k,1}';


        % Define var
        Q = zeros(1,(length(qmax)-1));
        qin = qearo{i,1};
        qout = qoutc{i,1};
        % Now determine the het load for all concepts, as well as the corrected
        % total heat flux.
        for j = 1:(length(qmax)-1)

            qtoto(j+1,1) = (qin(j)-qout(j));
            Q(j+1) = qtoto(j+1,1)*dt+Q(j);
            if Q(j+1) < 0
                Q(j+1) = 0;
                qtoto(j+1,1) = 0;
            end

        end
        qtot{i,1} = qtoto;
        Qtot{i,1} = Q;

        Qmax(i,1) = max(Q);
    end
    
    fraction(:,k) = Qmax./Qmax(4);
end
    % plotting thermal results
%     figure(2)
%     plot(t,qtoto,'r')
%     hold on
%     plot(t,qmax,'b')
%     figure(3)
%     plot(t,Q)
% 
%     Qmax
%     Qmax./Qmax(4)

% Temperature range vector
Twall = linspace(min,mout,(mout-min+1));

figure(1)
plot(Twall,fraction(1,:),'b')
hold on
grid on
plot(Twall,fraction(2,:),'-.b')
%plot(Twall,fraction(3,:),'r')
plot(Twall,fraction(4,:),'r')
plot(Twall,fraction(5,:),'-.r')
title('Heat load fractions')
xlabel('wall temperature [K]')
ylabel('Heat load fraction [-]')
legend('ballute','isotensoid','stacked torroid','tension cone')
axis([min mout 0.55 1.25  ])   
    
    
    
    
    

