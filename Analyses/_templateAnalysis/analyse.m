% ------------------------------------------------------------------------
% THE HILL - ANALYSIS TOOL
%
% HOW TO USE:
% 1. Make selection of mastersheet measurement id's in analysisSetup.xlsx
% 2. Run this file. (No editing necessary/recommended)
%
% OPERATION:
% Reads settings from "TheHill_analysisSetup.xlsx" in the current folder.
% Extracts data via the mastersheet ("TheHill_Mastersheet.xlsx"), returns
% data, shows and saves default figures.
%
% MvN 2020 - Dimple Aerospace BV
% ------------------------------------------------------------------------
clc; clear; close all;

%% Read settings from setup file

% Look for setup file
file = dir('*config*.xlsx');    % Allow some filename modification
if isempty(file)
    error('Could not find the setup file.');
end

% Load setup file
[~, ~, setup] = xlsread(file.name);

% Extract data from setup file
id  = cell2mat(setup(2:11,2));
idx = find(~isnan(id));
id  = id(idx);
N   = length(id);

% Check settings
if length(id) ~= length(unique(id))
    error('Duplicate measurement IDs specified.');
end

forceFileRead = zeros(N,1);
for i = 1:N
    spec = setup{1+i,3};
    if strcmp(spec,'yes')
        forceFileRead(i) = 1;
    end
end

colors  = setup(2:N+1,4);
lines   = setup(2:N+1,5);
markers = setup(2:N+1,6);

% Import data via mastersheet
for i = 1:N

    [name{i}, Cd0{i}, Re0{i}, dCd{i}, dCdp{i}, Re_target{i}, RMSE{i},...
        RMSE_X{i}, F{i}, F_rms{i}, F_power{i}, p{i}, T{i}, Troom{i},...
        pa{i}, hum{i}, Re{i}, V{i}, rho{i}, nu{i}, nu_avg{i}, corr{i},...
        y{i}, u{i}, u_rms{i}, u_power{i}, ut{i}, y0{i}, k{i}, B{i}, PI{i}, d{i},...
        d_star{i}, theta{i}, H{i}, up_model{i}, yp_model{i}] ...
        = importFromMastersheet(id(i), forceFileRead(i));
    
end

% Make default plots
defaultPlots(name, Cd0, Re0, dCd, dCdp, Re_target, RMSE, RMSE_X, F, F_rms,...
    F_power, p, T, Troom, pa, hum, Re, V, rho, nu, corr, colors, lines, markers);

BLplots(name, y, u, u_rms, u_power, nu_avg, ut, y0, k, B, PI, d, d_star,...
    theta, H, up_model, yp_model, V, dCdp);

% Save figures to file
% saveAllFigures

%% CUSTOM CODE

