close all;
x = linspace(0,1,21);
p = rand(5,1);
for i = 1:10
p = rand(10,1)-0.5;
% p = ones(3,1);
p(end-1:end) = 0;
hold on

plot(x,polyval(p,x)/sum(p)*1)
end
