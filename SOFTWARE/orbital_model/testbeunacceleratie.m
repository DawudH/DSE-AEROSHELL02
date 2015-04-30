V = 0;
a = 10;
s = 0;
dt = 0.1;
tend = 10;
t = 0:dt:tend;
s_real = 0.5*a*t.^2;
for tc = 1:length(t)-1
    V(tc+1) = V(tc) + a*dt;
    s(tc+1) = s(tc) + V(tc)*dt + 0.5*a*dt^2;
end

plot(t, s, 'green');
hold on;
plot(t, s_real);
s_real(end)-s(end)