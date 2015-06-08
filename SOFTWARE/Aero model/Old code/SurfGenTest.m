clear; close all
Nr = 40
for i =1:10
    h = rand(1);
    a = rand(1);
    r = 6;
    poly = -1 +2*rand(1,10);
    poly(end-1:end) = 0;
    TriGeom = ParaGeom(Nr,a,r,h,poly);
    subplot(5,2,i)
    trisurf(TriGeom)
end
axis equal
    
