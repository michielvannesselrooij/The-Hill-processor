function [T, Troom, pa, hum, tunnel_q, tunnel_T, tunnel_pa,...
    tunnel_hum, tunnel_rpm, p, q] = read_file_env(fileName, tunnel)
% ------------------------------------------------------------------------
% Reads various data from .csv measurement file
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    data = importdata(fileName);
    
    % Separate measurements 
    % A) via empty rows
    filter1 = [1; find(sum(abs(data),2) == 0)+1; length(data)+1];

    % B) via steps in rpm
    rpm = data(:,7);
    rpmChange = rpm(2:end) - rpm(1:end-1);
    filter2 = find(rpmChange > 5)+1;

    filter = sort(unique( [filter1; filter2] ));
    
    % Initialize
    T           = zeros(length(filter)-1,1);
    Troom       = T;
    pa          = T;
    hum         = T;
    tunnel_T    = T;
    tunnel_q    = T;
    tunnel_pa   = T;
    tunnel_hum  = T;
    tunnel_rpm  = T;
    p           = repmat(T,1);
    q           = T;
    
    % Extract
    j=1;
    for i = 1:length(filter)-1
        
        % Import Hill sensor data
        T(j)     = mean(data( filter(i) : filter(i+1)-1, 4));              % Temperature test section [C]
        Troom(j) = mean(data( filter(i) : filter(i+1)-1, 5));              % Temperature at control box [C]
        pa(j)    = mean(data( filter(i) : filter(i+1)-1, 6));              % Absolute pressure [Pa]
        hum(j)   = mean(data( filter(i) : filter(i+1)-1, 3));              % Relative humidity [%]
        q(j)     = mean(data( filter(i) : filter(i+1)-1, 7));              % Mensor (TEMP)
        
        
        % Optionally import tunnel sensor data
        if strcmp(tunnel,'M-tunnel')
        
            tunnel_q(j)   = mean(data( filter(i) : filter(i+1)-1, 7));     % Dynamic pressure (q)
            tunnel_T(j)   = mean(data( filter(i) : filter(i+1)-1, 8));     % Temperature test section [C]
            tunnel_pa(j)  = mean(data( filter(i) : filter(i+1)-1, 9))*100; % Absolute pressure [Pa]
            tunnel_rpm(j) = mean(data( filter(i) : filter(i+1)-1, 10));    % Fan rpm
            tunnel_hum    = hum(j);                                          % Not available. Use Hill data.
        
        % elseif strcmp(tunnel,'W-tunnel')
        
        % elseif strcmp(tunnel,'A-tunnel')
        
        % elseif strcmp(tunnel,'LTT')
        
        else
            
            tunnel_q   = [];
            tunnel_T   = [];
            tunnel_pa  = [];
            tunnel_rpm = [];
            tunnel_hum = [];
        
        end
        
        j=j+1;
        
    end

    disp(['Succesfully imported file ' fileName]);
    disp(' ');

end