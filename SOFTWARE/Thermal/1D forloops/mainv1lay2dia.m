%% Thermal tool FR
% Suthes & Lucas
% v1 - Implement contact resistance
clear all
close all
clc

%% Define input

% lay-up (t[mm], k[w/m/K], rho[kg/m3], cp[J/kg/K], emis[-], allowT[K])
filename   = 'layup2.txt';
layupin    = dlmread(filename);
% SI units
thicks = dlmread('thick.txt');
m = [5.0,0,0];
qfact = 1.0;
for c=10:1:15
    for d=160:1:180
        disp([num2str(c),', ',num2str(d)])
        layupin(1,1) = thicks(c);
        layupin(2,1) = thicks(d);
        layup      = zeros(size(layupin));
        layup(:,1) = layupin(:,1)./1000;
        layup(:,2:6) = layupin(:,2:6);
        L = layup(:,1);

        % Aero input, qsdot
        load('heatflux_out_d18_just_orbit','T','t','qmax_array')
        tq = find(not(qmax_array<0.01));
        qaero = qmax_array(tq(1):tq(end))*qfact;
        Tatm  = T(tq(1):tq(end)); 
        timeq = t(tq(1):tq(end))-t(tq(1));
        clear('T','t')

        % time and aeroheat
        ttot = timeq(end);  % end time of the orbit [s]
        dt   = 0.1;    % time step, chooseable
        T0 = 293;
        q0 = qaero*10000*qfact;
        nmax = int32(ttot/dt);  % number of time steps
        t = [0:double(nmax-1)]*dt;
        fact    = 1;           % multiplication factor of number of space steps
        kfact = [2.5e-5/fact;2.5e-6/fact;2.5e-2/fact]; % Conductivity factors


        % spacing
        L = int32(round(L*10000000));


        if length(L)==1
            gcd_old = L(1);
        else
            gcd_old = gcd(L(1),L(2));
            if length(L)>2
                for  j = 3:length(L)
                    gcd_old = gcd(gcd_old,L(j));
                end
            end
        end
        maxdx = double(gcd_old)/10000000;
        L = double(L)/10000000;


        Ltot = sum(L);      % length of the layup [m]
        dx   = maxdx/fact;     % length of a space step
        imax =  int32(Ltot/dx + 1); % maximum amount of points in spaceing

        % Radiation parameters
        emiss = layup(1,5);
        sig   = 5.670373e-8; %[W/m2/K4]

        %% Material properties
        layers = length(L)-1;
        k = zeros(imax+layers-1,1);
        rho = zeros(imax+layers-1,1);
        cp = zeros(imax+layers-1,1);
        x = zeros(imax+layers,1);

        %determine material layer thicknesses
        indexx = zeros(length(L)+1,1);
        indexx(1) = 1;
        indexx(2:end) = int32(cumsum(L/dx))+1;
        for j = 1:length(L);
            k(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,2);
            rho(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,3);
            cp(indexx(j)+j-1:indexx(j+1)+j-2) = layup(j,4);
            x(indexx(j)+j-1:indexx(j+1)+j-1) = dx*[indexx(j)-1:indexx(j+1)-1];
        end
        for j = 1:layers
            k(indexx(j+1)+j-1) = kfact(j)*(k(indexx(j+1)+j-2)+k(indexx(j+1)+j))/(k(indexx(j+1)+j-2)*k(indexx(j+1)+j));
            rho(indexx(j+1)+j-1) = (rho(indexx(j+1)+j-2)+rho(indexx(j+1)+j))/2;
            cp(indexx(j+1)+j-1) = (cp(indexx(j+1)+j-2)+cp(indexx(j+1)+j))/2;
        end
        alpha = k./rho./cp;
        v = alpha*dt/dx/dx;


        %% Implement Cranck-Nickelson with voids

        % i is space, n is time
        % space is rows, time is columns
        T = zeros(imax+layers,nmax);
        T(:,1) = T0;
        qs = interp1(timeq,qaero,t)*10000;
        Tamb = interp1(timeq,Tatm,t);

        % Define matrices for Crank-Nicolson
        w = zeros(imax+layers,1);
        w(1:end-1) = v/2;
        w(2:end) = w(2:end) + v/2;
        C = diag(1+w) - diag(0.5*v,1) - diag(0.5*v,-1);
        N = diag(1-w) + diag(0.5*v,1) + diag(0.5*v,-1);

        A = zeros(imax+layers,1); 
        Cinv = inv(sparse(C));
        G = full(Cinv*sparse(N));

        % Perform scheme over time
        for n=1:nmax-1
            qrs = -emiss*sig*(T(1,n)^4-Tamb(n)^4); % Radiation front
            qrb = -emiss*sig*(T(end,n)^4-T0^4); % Radiation back
            A(1) = ((qs(n)+qrs)/k(1))*v(1)*dx;
            A(end) = (qrb/k(end))*v(end)*dx;
            H = full(Cinv*sparse(A));
            T(:,n+1) = G*T(:,n) + H;
        end

        %% Output
        % Contour plots
        contourplot = 0;
        if contourplot
            xcont = x;
            for j = 1:length(L)-1
                xcont(indexx(j+1)+j-1:indexx(j+1)+j) = [x(indexx(j+1)+j-1)-dx/100,x(indexx(j+1)+j)+dx/100];
            end
            figure;
            hold on
            contourf(t,xcont,T)
            colormap parula
            colorbar
            title('Temperature over time and distance')
            xlabel('Time [s]')
            ylabel('Depth [m]')
            for j = 1:length(L)-1
                plot([0,nmax*dt],[x(indexx(j+1)+j),x(indexx(j+1)+j)],'--','Color', [255.0/256.0,130.0/256.0,28.0/256.0])
            end
            axis ij   
        end

        % Layer analysis
        results = zeros(length(L),1);
        layernames = cell(length(L),1);
        for j=1:length(L)
            results(j) = max(max(T(indexx(j)+j-1:indexx(j+1)+j-1,:)));
            layernames{j} = strcat('Layer',int2str(j));
        end
        output = table(layup(:,1)*1000,results,layup(:,6),layup(:,2),layup(:,3),layup(:,4),'RowNames',layernames,'VariableNames',{'Thickness','maxT','allowT','k','rho','cp'});

        if sum(results<layup(:,6))==4
            if sum(layup(:,1).*layup(:,3))<m(1)
                m = [sum(layup(:,1).*layup(:,3)),c,d];
            end
        end
    end
end
%disp(output)    
%disp(['m/A = ',num2str(sum(layup(:,1).*layup(:,3))),' [kg/m3]'])  
    
    
