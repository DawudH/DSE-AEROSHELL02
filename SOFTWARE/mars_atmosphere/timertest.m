clear all;
atm = marsatmosphere();
h = linspace(0,4e5,9000);
tic;
for i = h
    T = atm.getCheapDensity(i);
end
toc;