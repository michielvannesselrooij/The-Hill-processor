function [files] = getFileNames(folder)
% ------------------------------------------------------------------------
% Returns the master file names of all measurements in a specified folder
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    % Identify unique measurements
    files = dir([folder filesep '*_F.csv']);

    if ~isempty(files)
        % Get master file name
        for i=1:length(files)
            files2{i} = [files(i).name(1:end-6) '.csv'];
        end

    else
        % Get .mat files in case .csv files were moved/zipped
        files = dir([folder filesep '*.mat']);

        if ~isempty(files)
            % Get master file name (csv names expected by other functions)
            for i=1:length(files)
                files2{i} = [files(i).name(1:end-4) '.csv'];
            end
        
        % If those also don't exist, tell the user and abort if no hotwire    
        else
            if isempty(dir([folder filesep 'hotwire']))
                error(['There are no suitable drag measurement files in folder ' folder]);
                
            else
                warning(['There are no suitable drag measurement files in folder ' folder]);
            end
            
           files = [];
           return
           
        end

    end
    
    % Replace with final names
    files = files2;
    
    % Report measurements
    disp('The following measurements will be analyzed:');
    for i=1:length(files)
        disp(['  ' files{i}(1:end-4)]);
    end
    disp(' ');

end