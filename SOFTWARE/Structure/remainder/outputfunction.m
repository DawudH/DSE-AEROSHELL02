function [] = output(x1range,x2range,x1,x2,desired_output)

x1range = x1range;
x2range = x2range;

Cd,qmax,mf,M,T,...
Do,Di,theta,L,N,I,...
selected_material_gore,selected_material_rad,selected_material_axi,...
selected_material_gasb,selected_material_toro,typeoftoroid

a = 0.5;
b = 3000;
c = 2;
d = 22;
e = 273.15;
f= 12;
g= 5;
h = 60;
i = 1;
j = 8;
k = 18;

variables = { {'Cd'}, {'qmax'}, {'mf'}, {'M'}, {'T'}, {'Do'}, {'Di'},...
    {'theta'}, {'L'},{'N'},{'I'}};
values = [a;b;c;d;e;f;g;h;i;j;k]

idx = find(strcmp([variables{:}], x1));
a   = values(idx);



m_total_1 = zeros(length(x1range),length(x2range));
m_total_2 = zeros(length(x1range),length(x2range));
m_total_3 = zeros(length(x1range),length(x2range));

i = 1;
for vary_1=x1range
    j = 1;
    for vary_2=x2range
        
        idx_x2 = find(strcmp([variables{:}], x2));
        a = values(idx
        [mtoroidbar,mgasbar,mfactor,mtotalbar,mtotal,pmin] = ST_TC_TIAD_massfunction(a,b,c,d,e,f,g,h,i,j,k,...
        'Kevlar 49','Kevlar 49','Kevlar 49','Kevlar 49','Kevlar 49',0);
        
        m_total_1(i,j)=mtotal(1);
        m_total_2(i,j)=mtotal(2);
        m_total_3(i,j)=mtotal(3);
        j = j + 1;
    end
    
    i = i + 1;
end

Cdplot = repmat(x1range',[1 length(x2range)]);
Doplot = repmat(x2range',[1 length(x1range)]);

subplot(131)
mesh(x1range,x2range,m_total_1')
subplot(132)
mesh(x1range,x2range,m_total_2')
subplot(133)
mesh(x1range,x2range,m_total_3')