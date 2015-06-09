clc
clear
close all

addpath('..\..\Aero model\orbits')

load('out_d6_just_orbit','out')
outj6 = out;
load('out_d9_just_orbit','out')
outj9 = out;
load('out_d12_just_orbit','out')
outj12 = out;
load('out_d15_just_orbit','out')
outj15 = out;
load('out_d18_just_orbit','out')
outj18 = out;

load('out_d6_one_time','out')
outo6 = out;
load('out_d9_one_time','out')
outo9 = out;
load('out_d12_one_time','out')
outo12 = out;
load('out_d15_one_time','out')
outo15 = out;
load('out_d18_one_time','out')
outo18 = out;

% figure('name','radius')
% hold on
% just=fit([6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]','cubicinterp');
% plot(just,[6,9,12,15,18],[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)])
% one=fit([6,9,12,15,18]',[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)]','cubicinterp');
% plot(one,[6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)])
% 
% figure('name','radius')
% hold on
% just = polyfit([6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]',3);
% justx = linspace(6,18,100);
% justy = polyval(just,justx);
% plot(justx,justy,[6,9,12,15,18]',[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)]','x');
% one = polyfit([6,9,12,15,18]',[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)]',3);
% onex = linspace(6,18,100);
% oney = polyval(one,onex);
% plot(onex,oney,[6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)],'d')
cc = parula(7);

figure('name','radius')
hold on
plot([6,9,12,15,18],[max(outj6.q),max(outj9.q),max(outj12.q),max(outj15.q),max(outj18.q)],'-o','color',cc(1,:),'MarkerEdgeColor',cc(1,:))
plot([6,9,12,15,18],[max(outo6.q),max(outo9.q),max(outo12.q),max(outo15.q),max(outo18.q)],'-d','color',cc(3,:),'MarkerEdgeColor',cc(3,:))