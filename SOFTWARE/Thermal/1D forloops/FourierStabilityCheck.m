clear all
close all
clc

%% Determination of the amplification factor and Fourier stability analysis.

% Fourier wavenumvers
beta = 0:0.01:pi;

% Material properties Nextel BF-20
k   = 0.146;
rho = 1362;
cp  = 1130;

% Used coefficient
c = k/(rho*cp);

% Scheme refinement
CFL = [5e-5 1e-5];
dx  = 5e-4;
dt  = (CFL*dx)/(c); 

% Calculate the amplification factor
ampact = zeros(length(beta),length(CFL));
for i = 1:length(CFL)
    ampfact(:,i) = ((CFL/dx).*cos(beta)+(1-(CFL/dx))) ./ ((CFL/dx).*cos(beta)+(1+(CFL/dx)));
end

% Plotting
plot(beta,abs(ampfact(:,1)))
axis([0 pi 0.6 1.2])
grid on
title('Fourier amplification factor')
xlabel('beta')
ylabel('amplification')