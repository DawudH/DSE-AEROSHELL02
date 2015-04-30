begin = -4e6;
ending = -5e6;
step = -1e5;
refine = 10;
rx = begin:step:ending; %[m]
CD = 1.2; %[-]
fid = fopen('orbit.txt','a');
for CD=0.5:0.1:1.5
    find_orbit(rx,CD,fid,begin,ending,step,refine);
end
fclose(fid);