function [Cd, F, F_rms, F_power, p, T, Troom, pa, hum, Re, V, rho, nu,...
    nu_avg, corr, y, u, u_rms, u_power, ut, y0, k, B, PI, d, d_star, theta, H, ...
    up_model, yp_model] = importFolder(folder, pConfig, pCount, ...
    qChannel, forceFileRead, useTunnelData, Re_corr, warmUp)
% ------------------------------------------------------------------------
% Reads all measurements in the specified folder
% 
% INPUT:
% folder        : relative folder path
% forceFileRead : forces reading of .csv files even when .mat file exists
% pConfig       : pressure channel configuration id from excel file
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

% Check input
if ~exist('forceFileRead','var')
    forceFileRead = 0;
end

if ~exist('useTunnelData','var')
    useTunnelData = [];
end

% Find full path
current    = cd;                                             % Current path
baseFolder = [filesep 'Analyses'];                           % Analysis folder
idx        = strfind(current,baseFolder)-length(baseFolder); % Trim path idx
base       = current(1:idx+length(baseFolder)-1);            % Trim path

folder   = [base filesep 'Measurements' filesep folder];

% Retrieve file names    
files = getFileNames(folder);

% Import
F     = cell(size(files));
F_raw = F;
V     = F;

% Initialize force measurement outputs
Cd = cell(size(files));
F  = cell(size(files));
F_rms = cell(size(files));
F_power = cell(size(files));
p = cell(size(files));
T = cell(size(files));
Troom = cell(size(files));
pa = cell(size(files));
hum = cell(size(files));
Re = cell(size(files));
V = cell(size(files));
rho = cell(size(files));
nu = cell(size(files));
corr = cell(size(files));


% Import data
for i=1:length(files)

    fileName = [folder filesep files{i}];  

    [Cd{i}, F{i}, F_rms{i}, F_power{i}, p{i}, T{i}, Troom{i}, pa{i},...
        hum{i}, Re{i}, V{i}, rho{i}, nu{i}, corr{i}]...
        = importMeasurement(fileName, forceFileRead, pConfig,...
            pCount, qChannel, useTunnelData, Re_corr);

end

% Ignore first measurement if it is a warm up (unless file it is deleted)
if exist('warmUp','var')
    if warmUp == 1 && length(files)/2 == floor(length(files)/2)
        Cd(1)       = [];
        F(1)        = [];
        F_rms(1)    = [];
        F_power(1)  = [];
        p(1)        = [];
        T(1)        = [];
        Troom(1)    = [];
        pa(1)       = [];
        hum(1)      = [];
        Re(1)       = [];
        V(1)        = [];
        rho(1)      = [];
        nu(1)       = [];
        corr(1)     = [];
    end
end

% Import hotwire data if available
hotwireFolder = [folder filesep 'hotwire'];
HW_mat_files = dir([hotwireFolder filesep '*.mat']);
HW_csv_files = dir([hotwireFolder filesep '*HW.csv']);

if exist(hotwireFolder,'dir') && ~(isempty(HW_mat_files) && isempty(HW_csv_files))

    % Determine average viscosity
    if isempty(nu)
        
        nu_avg = 1.51e-5;
        warning('Using default value for nu');
        
    else
        
        nu_avg = 0;
        for i=1:length(nu)
            nu_avg = nu_avg+mean(nu{i});
        end
        nu_avg = nu_avg/length(nu);
        
    end
    
    
    if forceFileRead==0 && ~isempty(HW_mat_files)
        
        % load processed data if available
        for i=1:length(HW_mat_files)
            
            load([HW_mat_files(i).folder filesep HW_mat_files(i).name],...
                'y', 'u', 'u_rms', 'u_power', 'ut', 'y0', 'k', 'B', 'PI',...
                'd', 'd_star', 'theta', 'H', 'up_model', 'yp_model');
            
            y2{i}        = y;
            u2{i}        = u;
            u_rms2{i}    = u_rms;
            u_power2{i}  = u_power;
            ut2(i)       = ut;
            y02(i)       = y0;
            k2(i)        = k;
            B2(i)        = B;
            PI2(i)       = PI;
            d2(i)        = d;
            d_star2(i)   = d_star;
            theta2(i)    = theta;
            H2(i)        = H;
            up_model2{i} = up_model;
            yp_model2{i} = yp_model;
            
        end
        
        y        = y2;
        u        = u2;
        u_rms    = u_rms2;
        u_power  = u_power2;
        ut       = ut2;
        y0       = y02;
        k        = k2;
        B        = B2;
        PI       = PI2;
        d        = d2;
        d_star   = d_star2;
        theta    = theta2;
        H        = H2;
        up_model = up_model2;
        yp_model = yp_model2;
        
    elseif ~isempty(HW_csv_files)
        
        % Import data
        for i=1:length(HW_csv_files)
            
            disp('');
            disp(['Importing hotwire data: ' HW_csv_files(i).name]);
            disp('');
            
            dataFile   = [hotwireFolder filesep HW_csv_files(i).name];
            calFile    = [dataFile(1:end-4) '_cal.csv'];
            wallFile   = dir([dataFile(1:end-18) '*wall.csv']);
            ignoreFile = dir([dataFile(1:end-18) '*ignore*']);
            
            % Use the first calibration file in folder, if it doesn't exist
            % for this case
            if ~exist(calFile, 'file')
                calibrationFiles = dir([hotwireFolder filesep '*HW_cal.csv']);
                if length(calibrationFiles)<1
                    % No calibration files in folder
                    error(['No calibration files for the hotwire'...
                        'measurements in ' hotwireFolder]);
                else
                    % Otherwise use the first calibration file in folder
                    calFile  = [calibrationFiles(1).folder filesep...
                                        calibrationFiles(1).name];
                end
            end 
            
            % Find a wall calibration if it doesn't exist for this case
            if ~exist(wallFile, 'file')
                wallCalFiles = dir([hotwireFolder filesep '*wall.csv']);
                if length(wallCalFiles)<1
                    % No calibration files in folder
                    % warning(['No wall calibration files for the hotwire'...
                    %    'measurements in ' hotwireFolder]);
                else
                    % Otherwise use the first calibration file in folder
                    wallFile  = [wallCalFiles(1).folder filesep...
                                        wallCalFiles(1).name];
                    warning(['Wall calibration ' wallCalFiles(1).name ...
                        ' used for ' HW_csv_files(i).name]);
                end
                
            else
                wallFile = [wallFile.folder filesep wallFile.name];
                
            end

            % Determine data points to ignore based on ignore .txt file
            ignore = [];
            if ~isempty(ignoreFile)
                str = erase(ignoreFile.name, '.txt');
                idx = strfind(str, 'ignore_') + 7;
                ignore = str2double(str(idx:end));
            end
            
            % Process data
            [y{i}, u{i}, u_rms{i}, u_power{i}, ut(i), y0(i), k(i), B(i),...
                PI(i), d(i), d_star(i), theta(i), H(i), up_model{i}, ...
                yp_model{i}] = hotwireProcessor(dataFile, ...
                calFile, wallFile, ignore, nu_avg);
            
        end
    end
    
else
    
    % Leave outputs empty if no hotwire data is found
    y        = [];
    u        = [];
    u_rms    = [];
    u_power  = [];
    ut       = [];
    nu_avg   = [];
    y0       = [];
    k        = [];
    B        = [];
    PI       = [];
    d        = [];
    d_star   = [];
    theta    = [];
    H        = [];
    up_model = [];
    yp_model = [];
    
end