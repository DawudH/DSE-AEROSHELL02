clc
clear all
close all

% load constants
constants

ry = 10*R_m; %[m]
v = 7000; %[m/s]


dt_atmos = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8];
dt_init = dt_atmos;
dt_kep_init = 1e-8;
rx = 4190000;
tend = 3600 * 24* 1;

control.CL_range = [-0.35 0.35];
control.CLCD = 0.25;
control.a = 2.9*g_earth;
control.CLa = 0.02;
control.dalpha = 0.2;
control.CL_init = -0.18;

% run for different timesteps
R = cell(length(dt_atmos),1);
error = cell(length(dt_atmos),1);
max_error = zeros(length(dt_atmos),1);
    for i = 1:length(dt_atmos)
        [out] = orbit_full(rx,ry,v,dt_init(i),dt_atmos(i),dt_kep_init,tend,control);
        R{i} = out.R;
        len_old = length(R{1});
        len_new = length(out.R);
        if i > 1
            k = 1:2^(i-1):len_old;
            error{i} = out.R(1:length(k),:) - R{1}(k,:);
        else
            error{i} = out.R - out.R;
        end
        max_error(i,1) = max(sqrt(error{i}(:,1).^2 + error{i}(:,2).^2 + error{i}(:,3).^2));
        
    end
    
    figure('name','discretization error')
    grid on
    plot(dt_atmos,max_error)
    
    
