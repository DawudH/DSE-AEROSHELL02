clear

n = 10;
vali = ones(n,n,n);
vali(1,1,1) = 0;
go = 0;
dx = linspace(-5,5,n);
dz = linspace(-5,5,n);
dy = linspace(-5,5,n);

x = 300;
y = 1;
z = 100;

Mx = 0;
My = 0;
Mz = 0;

Mxcg = zeros(n,n,n);
Mycg = zeros(n,n,n);
Mzcg = zeros(n,n,n);
while (max(max(max(vali)))-min(min(min(vali))))>0.001
    for i = 1:length(dx)
        for j = 1:length(dy)
             for k = 1:length(dz)
                Mxcg(i,j,k) = Mx + y*dz(k) - z*dy(j);
                Mycg(i,j,k) = My - x*dz(k) + z*dx(i);
                Mzcg(i,j,k) = Mz - y*dx(i) + x*dy(j);
             end
        end
    end
    Mcg = abs(Mxcg) + abs(Mycg) + abs(Mzcg);
    [vali,rowi] = min(Mcg);
    [valj,rowj] = min(vali);
    [valk,rowk] = min(valj);
     loc = [rowi(1,rowj(1,1,rowk),rowk),rowj(1,1,rowk),rowk];
     if max(loc)<n && min(loc)>1
         dx = linspace(dx(loc(1)-1),dx(loc(1)+1),n);
         dy = linspace(dy(loc(2)-1),dy(loc(2)+1),n);
         dz = linspace(dz(loc(3)-1),dz(loc(3)+1),n);
     elseif min(loc)==1
         dx = linspace(dx(loc(1)),dx(loc(1)+1),n);
         dy = linspace(dy(loc(2)),dy(loc(2)+1),n);
         dz = linspace(dz(loc(3)),dz(loc(3)+1),n);
     elseif max(loc)==n
         dx = linspace(dx(loc(1)-1),dx(loc(1)),n);
         dy = linspace(dy(loc(2)-1),dy(loc(2)),n);
         dz = linspace(dz(loc(3)-1),dz(loc(3)),n);
     end
     go = go+1
end
 dxv = dx(n/2)
 dyv = dy(n/2)
 dzv = dz(n/2)
