function [p, p_rms, p_power, p_raw] = read_file_p(fileName, sampleTime)
% ------------------------------------------------------------------------
% Reads Evoscann pressure data from .csv measurement file
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------
        
    data   = importdata(fileName);                      % import
    %data   = data*100;                                  % mbar to Pa
    
    % Clean NaN if present
    data(isnan(data)) = 0;

    % Separate measurements via empty rows
    filter = [0; find(sum(abs(data),2) == 0)];
    p_raw  = cell((length(filter)-1),size(data,2));     % initialize
    p      = zeros(size(p_raw));                        % initialize

    for i = 1:size(p,1)
        for j = 1:size(p,2)
            
            p_raw{i,j}   = data( filter(i)+1 : filter(i+1)-1, j);
            
            % Exclude missing values
            p_raw{i,j}(p_raw{i,j}==0) = [];
            
            p(i,j)       = mean(p_raw{i,j});
            p_rms(i,j)   = rms(p_raw{i,j});
            
            % Power spectrum
            n        = length(p_raw{i,j});  % Number of samples
            fs       = n/sampleTime;        % Sample rate

            y        = fft(p_raw{i,j}-mean(p_raw{i,j}));
            f        = (0:n-1)*(fs/n);      % Frequency range vector
            power    = abs(y).^2/n;         % Power
            
            f        = f(1:floor(n/2));     % Second half is symmetric
            power    = power(1:floor(n/2)); % Second half is symmetric
    
            if sum(power>0) > 0
                [power, f] = findpeaks(power,f','MinPeakProminence',0.1);
                p_power{i,j} = [f, power];      % Store in output
            else
                p_power{i,j} = [];
            end

        end
    end
    
    if isempty(p)
        p       = 0;
        p_rms   = 0;
        p_power = 0;
        p_raw   = 0;
    end

    % Tell what I've done
    disp(['Succesfully imported file ' fileName]);
    disp(' ');

end