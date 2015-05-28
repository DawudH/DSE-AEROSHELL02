clear all
close all
x = [0 -0.5 -0.25 0 0.25 0.5 0]; y = [0 0 0.5 1 0.5 0 0];
s45 = 1/sqrt(2); 
w =[1 0.5 1 0.5 1 0.5 1];
circle = rsmak([0 0 0 1/3 1/3 2/3 2/3 1 1 1], [w.*x;w.*y;w]);

fnplt(circle)
axis equal
hold on
plot(circle.coefs(1,:),circle.coefs(2,:),'x')
plot(x,y,'o')
