dt = [0.01,0.1,0.2,0.5,1];
e = [5.5, 6.5, 7.3, 10.0, 13.7];

figure('name','Error between Kepler and numerical simulation');
loglog(dt,e,'d-');
grid on