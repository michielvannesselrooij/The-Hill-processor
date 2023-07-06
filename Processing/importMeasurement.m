function [Cd, F, F_rms, F_power, p, P, X, T, Troom, pa, hum, Re, V, rho, nu, corr] = ...
    importMeasurement(fileName, forceFileRead, pConfig, useTunnelData, Re_corr)
% ------------------------------------------------------------------------
% Reads all data from .csv measurement files and stores it as .mat, or 
% reads from .mat file if available and not overridden by user.
% 
% INPUT:
% fileName      : measurement csv file path name (exclude '_F' or '_p')
% forceFileRead : forces reading of .csv files even when .mat file exists
% pConfig       : pressure channel configuration id from excel file
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Import -----------------------------------------------------------------
% ------------------------------------------------------------------------

% Check if .mat file exists already
fileNameMat = [fileName(1:end-4) '.mat'];
previouslyProcessed = exist(fileNameMat,'file');

% Determine whether to import from .csv (slow) or .mat (fast)
if forceFileRead || previouslyProcessed~=2

    fileName_F = [fileName(1:end-4) '_F.csv'];
    fileName_p = [fileName(1:end-4) '_p.csv'];

    % Load data from .csv files
    [F, F_rms, F_power, sampleTime] = read_file_F(fileName_F);
    [T, Troom, pa, hum, tunnel_q, tunnel_T, tunnel_pa, tunnel_hum, ...
        tunnel_rpm, p, q_Mensor] = read_file_env(fileName, useTunnelData);

    if length(F) ~= length(T)
        error('The number of condition measurements does not match the number of force measurements');
    end

    if length(F) ~= size(p,1)
        error('The number of pressure measurements does noet match the number of force measurements');
    end
    
    % Calibrate force data
    showPlot = 1;
    F = sensorCalibration(F, showPlot);
    
    % Read pressures
    if exist(fileName_p, 'file')
       
        % Read pressure data file
        [p, p_rms, p_power] = read_file_p(fileName_p, sampleTime);
       
        % Interpret pressure data
        [q, dF_gap, dF_int, P, X, p_gap] = pressureTranslator(p, pConfig);
       
    end

    % Calculate Reynolds etc
    if useTunnelData
        [Re, V, rho, nu] = velocityCalculator(tunnel_q, tunnel_pa, tunnel_T, tunnel_hum);
    else
        [Re, V, rho, nu] = velocityCalculator(q, pa, T, hum);
    end
    
    % Save as .mat file
    save(fileNameMat, 'F', 'F_rms', 'F_power', 'p', 'p_rms', 'p_power',...
        'q', 'q_Mensor', 'T', 'Troom', 'pa', 'hum', 'tunnel_q', 'tunnel_T',...
        'tunnel_pa', 'tunnel_hum', 'tunnel_rpm', 'Re', 'V', 'rho', 'nu');

else
    
    % Load data from .mat file
    load(fileNameMat, 'F', 'F_rms', 'F_power', 'p', 'p_rms', 'p_power',...
        'q', 'q_Mensor', 'T', 'Troom', 'pa', 'hum', 'tunnel_q', 'tunnel_T',...
        'tunnel_pa', 'tunnel_hum', 'tunnel_rpm', 'Re', 'V', 'rho', 'nu');
    
    % Substitute in missing data if necessary
    if ~exist('F_rms','var')
        F_rms = 0;
    end
    if ~exist('F_power','var')
        F_power = 0;
    end
    
    % Interpret pressure data
    [~, dF_gap, dF_int, P, X, p_gap] = pressureTranslator(p, pConfig);
    
end

% Calculate pressure contribution to force
F_p = dF_gap + dF_int;

% Return Hill sensor data or tunnel sensor data?
if useTunnelData
    q   = tunnel_q;
    T   = tunnel_T;
    pa  = tunnel_pa;
    hum = tunnel_hum;
end

% --------------------------------------------------------------------
% Perform data corrections -------------------------------------------
% --------------------------------------------------------------------

% Correct q based on pitot position
q      = q - q(end)*q./max(q);
q      = q*(1+Re_corr);

% Null shift correction of force measurement
shift  = -(F(end)-F(1)) * (F-F(1))./(max(F)-F(1));
F      = F + shift;
F(end) = F(1);

% Combine corrections to pass to function output
corr = {shift, dF_gap, p_gap, dF_int};

% --------------------------------------------------------------------
% Final data preparation ---------------------------------------------
% --------------------------------------------------------------------

% Rerun velocity calculation
[Re, V, rho, nu] = velocityCalculator(q, pa, T, hum);

% Calculate C_D
S             = 0.3663 * 0.8813;    % test plate surface area
q             = 0.5*rho.*V.^2;

Cd.F          =  (F-F(1)) ./ (q.*S);
Cd.p          = -(F_p-F_p(1)) ./ (q.*S);
Cd.total      = Cd.F + Cd.p;

Cd.total(1)   = 0;                  % replace first value which is always NaN
Cd.F(1)       = 0;                  % replace first value which is always NaN
Cd.p(1)       = 0;                  % replace first value which is always NaN

% Replace values in case of missing data (e.g. pa/T/hum is zero)
Cd.total(isnan(Cd.total)) = 0;
Cd.F(isnan(Cd.F))         = 0;
Cd.p(isnan(Cd.p))         = 0;

Cd.total(isinf(Cd.total)) = 0;
Cd.F(isinf(Cd.F))         = 0;
Cd.p(isinf(Cd.p))         = 0;

Re(isnan(Re))             = 0;
Re(isinf(Re))             = 0;

% --------------------------------------------------------------------
% Error checks -------------------------------------------------------
% --------------------------------------------------------------------

N = length(F);
if size(p,1) ~= N
    warning(['Length of pressure data (' num2str(size(p,1)) ') does ' ...
        'not match length of F (' num2str(size(F,1)) ')']);
elseif size(pa) ~= N
    warning(['Length of environmental data (' num2str(size(pa,1)) ')' ...
        ' does not match length of F (' num2str(size(F,1)) ')']);
end