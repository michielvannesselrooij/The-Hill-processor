function [y, u, u_rms, u_power, ut, y0, k, B, PI, d, d_star, theta, H,...
    up_model, yp_model] = hotwireProcessor(dataFile, calFile, wallFile,...
    ignore, nu)

% -------------------------------------------------------------------------
% Hotwire settings
% -------------------------------------------------------------------------

% Probe temperature
Tw    = 225;        

% -------------------------------------------------------------------------
% Import calibration file
% -------------------------------------------------------------------------

% Import
delimiter = '\t';
startRow = 1;
formatSpec = '%s%f%[^\n\r]';
fileID = fopen(calFile,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '',...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
data = dataArray{2};

% Extract
filter0 = find(isnan(data));
filter1 = filter0(1:2:end-1);
filter2 = filter0(2:2:end);
u_cal   = str2double(dataArray{1}(filter1)); % Pitot velocities
T_cal   = str2double(dataArray{1}(filter2)); % Temperatures

volt_cal = zeros(length(filter1),1); 
filter1  = [filter1; length(data)+1];
for i = 1:length(filter1)-1
    raw         = data(filter1(i)+3 : filter1(i+1)-1);
    filtered    = raw( raw > mean(raw)-3*std(raw) & ...
                       raw < mean(raw)+3*std(raw));
    volt_cal(i) = mean(filtered);
end

% Calibration function (based on Hultmark and Smits, 2010)
% u_cal(isnan(u_cal))  = 0;
% [volt_cal, idx]      = sort(volt_cal);
% u_cal                = u_cal(idx);
% volt_cal(u_cal==Inf) = [];              % Catch potential corrupt data
% u_cal(u_cal==Inf)    = [];              % Catch potential corrupt data
% 
% T0     = 273.15;
% nu_cal = 1./pa_cal *4.18528*1e-4 .*(T_cal+T0).^2.5 ./(110.4+T_cal+T0);
% k_cal  = 418.4*( 5.75e-5*(1 + 0.00317*T_cal - 0.0000021*T_cal.^2) );
% y_cal  = volt_cal.^2 ./(k_cal.*(Tw-T_cal));
% 
% [p, ~, mu] = polyfit(y_cal, u_cal./nu_cal, 4);

% Create calibration 4-th order polynomial (based on convection theory)
u_cal(isnan(u_cal))  = 0;
[volt_cal, idx]      = sort(volt_cal);
u_cal                = u_cal(idx);
volt_cal(u_cal==Inf) = [];              % Catch potential corrupt data
u_cal(u_cal==Inf)    = [];              % Catch potential corrupt data

[p, ~, mu] = polyfit(volt_cal, u_cal, 4);

% -------------------------------------------------------------------------
% Import wall calibration file
% -------------------------------------------------------------------------

% n       = 4;            % Order of calibration fit
% p_wall  = zeros(1,n+1);
% mu_wall = zeros(2,1);

if ~isempty(wallFile)
    
    % Import
    fileID = fopen(wallFile,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '',...
        'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
        'TextType', 'string', 'ReturnOnError', false);
    fclose(fileID);
    data = dataArray{2};
    
    % Filter
    filter0 = find(isnan(data));
    filter1 = filter0(1:3:end-2);
    filter3 = filter0(3:3:end);
    
    % Extract
    T_wall  = str2double(dataArray{1}(filter3));        % Temperatures
    y_wall  = str2double(dataArray{1}(filter1))/1000;   % y values
    
    volt_wall = zeros(length(filter1),1); 
    filter1   = [filter1; length(data)+1];
    
    for i = 1:length(filter1)-1
        raw          = data(filter1(i)+3 : filter1(i+1)-1);
        filtered     = raw( raw > mean(raw)-3*std(raw) & ...
                            raw < mean(raw)+3*std(raw));
        volt_wall(i) = mean(filtered);
    end
%     volt_wall = volt_wall - mean(volt_wall(1:5));
    volt_wall = volt_wall - volt_wall(1);
    
    % Correct for temperature (Bruun, 1995)
    Tr = interp1(volt_cal, T_cal, volt_wall, 'linear', 'extrap');
    volt_wall = volt_wall.*sqrt((Tw-Tr)./(Tw-T_wall));
    
    % Curve-fit the wall correction 
    wallCorr = hotwireWallCalibration(volt_wall, y_wall);
    
    % [p_wall, ~, mu_wall] = polyfit(y_wall, volt_wall, n);
    
end

% -------------------------------------------------------------------------
% Import measurement file
% -------------------------------------------------------------------------

fileID = fopen(dataFile,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '',...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
data = dataArray{2};

% Extract
filter0 = find(isnan(data));
filter1 = filter0(1:3:end-2);
filter2 = filter0(2:3:end-1);
filter3 = filter0(3:3:end);
y       = str2double(dataArray{1}(filter1))/1000; % y values
u_inf   = str2double(dataArray{1}(filter2));      % Pitot velocities
T       = str2double(dataArray{1}(filter3));      % Temperatures
volt    = cell(length(filter1),1);                       
u_raw   = cell(length(filter1),1);  
u_power = cell(length(filter1),1);  
u       = zeros(length(filter1),1); 
u_rms   = zeros(length(filter1),1);  

% Calculate wall correction
if exist('volt_wall','var')
    for i=1:length(y)
        if y(i)<= max(y_wall)
            %wallCorr(i) = polyval(p_wall, y(i), [], mu_wall);
            wallCorr(i) = interp1(y_wall, volt_wall, y(i),...
                'linear', 'extrap');
        else
            wallCorr(i) = 0;
        end
    end
       
    figure; hold on; box on;
    [~, name, ~] = fileparts(dataFile);
    title(['Hot-wire wall correction for ' name], 'Interpreter', 'none');
    plot(volt_wall, y_wall, 'k-')
    plot(wallCorr, y, 'ks')
    
else 
    wallCorr = zeros(size(y));
end
    
% nu_HW   = 1./pa_cal *4.18528*1e-4 .*(T+T0).^2.5 ./(110.4+T+T0);
% k_HW    = 418.4*( 5.75e-5*(1 + 0.00317*T - 0.0000021*T.^2) );
    
filter1  = [filter1; length(data)+1];
for i = 1:length(filter1)-1
    
    % Raw voltage
    raw        = data(filter1(i)+3 : filter1(i+1)-1);
    volt{i}    = raw( raw > mean(raw)-3*std(raw) & ...
                      raw < mean(raw)+3*std(raw));
    
    % Correct for temperature (Bruun, 1995)
    Tr = interp1(volt_cal, T_cal, mean(volt{i}), 'linear', 'extrap');
    Tc = T(i);
    volt{i} = volt{i}*sqrt((Tw-Tr)/(Tw-Tc));

    % Apply wall correction
    volt{i} = volt{i} - wallCorr(i);

    % Translate voltage to velocity
%     u_raw{i}   = nu_HW(i)*polyval(p, volt{i}.^2/(k_HW(i)*(Tw-T(i))), [], mu);
    u_raw{i}   = polyval(p, volt{i}, [], mu);   % Raw velocity
    
    % Correct for freestream velocity error
    u_raw{i}   = u_raw{i} * u_inf(1)/u_inf(i);
    
    u(i)       = mean(u_raw{i});                % Mean velocity
    u_rms(i)   = rms(u_raw{i}-u(i));            % RMS velocity

    % Power spectrum
    n        = length(u_raw{i});               % Number of samples
    fs       = 10000;                          % Sample rate
    yf       = fft(u_raw{i}-mean(u_raw{i}));   % Fourier transform
    f        = (0:n-1)*(fs/n);                 % Frequency range vector
    power    = abs(yf).^2/n;                   % Power
    f        = f(1:floor(n/2));                % Second half is symmetric
    power    = power(1:floor(n/2));            % Second half is symmetric

    u_power{i} = [f', power];
    
end

u       = flipud(u);
y       = flipud(y);
u_raw   = flipud(u_raw);
u_rms   = flipud(u_rms);
u_power = flipud(u_power);


% -------------------------------------------------------------------------
% Remove lowest data points if necessary
% -------------------------------------------------------------------------
if ~isempty(ignore)
    u(1:ignore)         = [];
    y(end-ignore+1:end) = [];
    u_raw(1:ignore)     = [];
    u_rms(1:ignore)     = [];
    u_power(1:ignore)   = [];
end
% -------------------------------------------------------------------------
% Profile characterization
% -------------------------------------------------------------------------
showPlot = 0;
[ut, y0, k, B, PI, d, yp_model, up_model] = canonicalOpt(y, u, nu, showPlot);

U = u(end);
[~, d_star, theta, H] = calcBLproperties(u, y, U, d);


% -------------------------------------------------------------------------
% Store data
% -------------------------------------------------------------------------
save([dataFile(1:end-4) '.mat'], 'y', 'u', 'u_rms', 'u_power', 'ut', ...
    'y0', 'k', 'B', 'PI', 'd', 'd_star', 'theta', 'H', 'up_model', 'yp_model');