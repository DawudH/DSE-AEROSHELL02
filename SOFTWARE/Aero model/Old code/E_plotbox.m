function [] = E_plotbox( y, boxes, halfchord, values)
% [] = plotbox( y, boxes, halfchord, values)
%   3D-plot the wing box



boxes = boxes(1:10:end,:,1:10:end);
halfchord = halfchord(1:10:end);
values = values(1:10:end,1:10:end);


n = round(length(y)/10);
om = round(length(boxes(:,1))/10);
% size(boxes)
% size(halfchord)
% size(values)

for i = 1:n
    boxes(:,1,i) = boxes(:,1,i) + halfchord(i);
end


valuessize = size(values);

if valuessize(2) == 1
    boxessize = size(boxes);
    new = zeros(boxessize(3), boxessize(1));
    for i = 1:length(values);
        new(i,:) = values(i)*ones(boxessize(1),1);
    end
    values = new;
end


X = zeros(4, om*(n-1));
Y = X;
Z = X;
C = X;

for i = 1 : n-1
    for j = 1 : om
        Y(:,(om)*(i-1) + j) = [y(i), y(i), y(i+1), y(i+1)];
        if j < om
            X(:,(om)*(i-1) + j)=  [boxes(j,1,i), boxes(j+1,1,i), boxes(j+1,1,i+1), boxes(j,1,i+1)];
            Z(:,(om)*(i-1) + j)=  [boxes(j,2,i), boxes(j+1,2,i), boxes(j+1,2,i+1), boxes(j,2,i+1)];
            C(:,(om)*(i-1) + j)=  [values(i,j), values(i,j+1), values(i+1,j+1), values(i+1,j)];
        else
            X(:,(om)*(i-1) + j)=  [boxes(j,1,i), boxes(1,1,i), boxes(1,1,i+1), boxes(j,1,i+1)];
            Z(:,(om)*(i-1) + j)=  [boxes(j,2,i), boxes(1,2,i), boxes(1,2,i+1), boxes(j,2,i+1)];
            C(:,(om)*(i-1) + j)=  [values(i,j), values(i,1), values(i+1,1), values(i+1,j)];
        end
    end
end
f = fill3(X, Y, Z, C);
%set(f,'edgecolor', 'none');
minimum = min([min(min(X)), min(y), min(min(Z))]);
maximum = max([max(max(X)), max(y), max(max(Z))]);
%axis([minimum, maximum, minimum, maximum, minimum, maximum, min(min(C)), max(max(C))]);

end