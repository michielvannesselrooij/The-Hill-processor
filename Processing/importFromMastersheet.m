function [name, Cd0, Re0, dCd, dCdp, Re_target, RMSE, RMSE_X, F, F_rms,...
    F_power, p, T, Troom, pa, hum, Re, V, rho, nu, nu_avg, corr, y, u, ...
    u_rms, u_power, ut, y0, k, B, PI, d, d_star, theta, H, up_model, yp_model] ...
    = importFromMastersheet(id, forceFileRead)
% ------------------------------------------------------------------------
% Read the set of measurements referenced by the specified id in the
% measurement mastersheet. See TheHill_Mastersheet.xlsx.
% 
% INPUT:
% id            : identifier (number) found in the mastersheet
% forceFileRead : forces reading of .csv files even when .mat files exist
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

% Find mastersheet location
current    = cd;                                             % Current path
baseFolder = '/Analyses';                                    % Analysis folder
idx        = strfind(current,baseFolder)-length(baseFolder); % Trim path idx
base       = current(1:idx+length(baseFolder)-1);            % Trim path

file   = [base '/Measurement files/TheHill_Mastersheet.xlsx'];
        
% Import mastersheet data
[~, ~, master] = xlsread(file, 'Measurements');
        
% Select measurement
headers = 3;
idx     = find(cell2mat(master(headers+1:end,1))==id);
data    = master(headers+idx,:);

% Prepare settings
name            = data{3};  % (info) Name identifier
refPlate        = data{4};  % (info) ID of reference test plate
targetPlate     = data{5};  % (info) ID of target test plate
warmUp          = data{6};  % Disregard a first warm-up measurement?
folder          = data{8};  % Folder within default measurements folder
tunnel          = data{9};  % (info) What tunnel was the test performed in?
tunnelData      = data{10}; % Is tunnel data present in env file?
useTunnelData   = data{11}; % Should this data be used for velocity?
pConfig         = data{12}; % Pressure channel layout configuration
pCount          = data{13}; % Number of pressure channels
qChannel        = data{14}; % Pitot pressure channel
Re_corr         = data{15}; % Reynolds correction factor
if isnan(Re_corr)
    Re_corr = 0;
end

if ~exist('forceFileRead','var')
    forceFileRead = 0;
end

if strcmp(useTunnelData,'yes')
    useTunnelData = tunnel; % Use tunnel name to indicate settings
else
    useTunnelData = [];
end

if strcmp(warmUp,'yes')
    warmUp = 1;
else
    warmUp = 0;
end

% Import the folder with measurements
[Cd, F, F_rms, F_power, p, T, Troom, pa, hum, Re, V, rho, nu,...
    nu_avg, corr, y, u, u_rms, u_power, ut, y0, k, B, PI, d, d_star, theta, H, ...
    up_model, yp_model] = importFolder(folder, pConfig, pCount, ...
    qChannel, forceFileRead, useTunnelData, Re_corr, warmUp);

% Calculate delta's
if length(Cd) > 1
    
    [Cd0, Re0, dCd, dCdp, Re_target] = dragDelta(Cd, Re);

    % Calculate spread (RMSE)
    [RMSE, RMSE_X] = calcRMSE(dCdp, Re_target);
    
else
    Cd0 = []; Re0 = []; dCd = []; dCdp = []; Re_target = [];
    RMSE = []; RMSE_X = [];
end
