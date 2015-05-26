clear all
close all
x = [1  0 -1 0 1]; y = [0  1 0 -1 0];
s45 = 1/sqrt(2); 
w =[1 1 1 1 1];
circle = rsmak(augknt(0:3,3,1), [w.*x;w.*y;w])

fnplt(circle)
axis equal
hold on
plot(circle.coefs(1,:),circle.coefs(2,:),'x')
plot(x,y,'o')
