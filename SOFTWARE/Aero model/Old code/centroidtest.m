clear ; close all
for i = 1:10
x = linspace(-i,i,21);
y = (x.^2)/(i^2);
hold on
plot(x,y)
x1 = [0 x];
x1 = x1(1:end-1);
deltax = x-x1;

y1 = [0 y];
y1 = y1(1:end-1);
deltay = y-y1;

xover2 = x(1:end-1)+ deltax(2:end)/2;
yover2 = y(1:end-1)+ deltay(2:end)/2;

elelength = sqrt(deltax(2:end).^2+deltay(2:end).^2);

y_c(i) = sum(yover2.*elelength)/sum(elelength);
x_c(i) = sum(xover2.*elelength)/sum(elelength);

% hold on 
% plot(x,y)
% plot(xover2,yover2,'o')
% plot(x_c,y_c,'x')
% axis equal
end

