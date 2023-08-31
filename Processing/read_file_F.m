function [F, F_rms, F_power, sampleTime] = read_file_F(fileName)
% ------------------------------------------------------------------------
% Reads force data from .csv measurement file
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

startFileLoading = tic;

delimiter = '\t';
startRow = 4;
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(fileName,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);

% Measurement separating lines
filter = [0; find(dataArray{2}=='Y'); length(dataArray{2})+3];

F_raw = cell(length(filter)-1,1);   % initialize
F     = zeros(length(filter)-1,1);  % initialize

for i = 1:length(filter)-1

    items = filter(i)+1 : filter(i+1)-3;
    F_raw{i} = zeros(size(items));  % initialize

    k = 1;
    for j = items

        F_raw{i}(k) = str2double(dataArray{2}(j));
        k = k+1;

    end

    % Find mean (correct for truncation if needed)
    tolerance = 1e-4;
    timeout = 5;
    showPlot = false;
    F(i) = MeanWithBottomTruncation(F_raw{i}, tolerance, timeout, showPlot);
    
    % Power spectrum
    F_rms(i) = rms(F(i));

    fs           = 25000;               % Sample rate
    n            = length(F_raw{i});    % Number of samples
    
    y            = fft(F_raw{i}-mean(F_raw{i}));
    f            = (0:n-1)*(fs/n);     % Frequency range vector
    power        = abs(y).^2/n;         % Power        
    
    f            = f(1:floor(n/2));     % Second half is symmetric
    power        = power(1:floor(n/2)); % Second half is symmetric
    [pks, f_pks] = findpeaks(power,f,'MinPeakProminence',0.1); % Peaks
    
    idx          = find(f<100);         % Keep only data below 100Hz
    f            = f(idx);              % Second half is symmetric
    power        = power(idx);          % Second half is symmetric
    
    F_power{i}   = {[f', power'], [f_pks', pks']}; % Store in output
     
    % Deduce measurement time (e.g. for pressure data interpretation)
    sampleTime = n/fs;

end

% Tell what I've done
disp(['Succesfully imported file ' fileName]);
toc(startFileLoading);
disp(' ');
