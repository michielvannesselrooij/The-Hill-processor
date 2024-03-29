function [T, Troom, pa, hum, tunnel_q, tunnel_T, tunnel_pa,...
    tunnel_hum, tunnel_rpm, p] = read_file_env(fileName, tunnel)
% ------------------------------------------------------------------------
% Reads various data from .csv measurement file
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    data = readmatrix(fileName);
    
    % Exclude any empty rows (legacy data)
    zeroRows = find(sum(abs(data),2) == 0);
    data(zeroRows, :) = [];

    % Separate measurements by NaN in temperature
    filter = [1; find(isnan(data(:,2))); size(data,1)+1];
    
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
        T(j)     = mean(data( filter(i)+1 : filter(i+1)-1, 1));              % Temperature test section [C]
        Troom(j) = mean(data( filter(i)+1 : filter(i+1)-1, 2));              % Temperature at control box [C]
        pa(j)    = mean(data( filter(i)+1 : filter(i+1)-1, 4));              % Absolute pressure [Pa]
        hum(j)   = mean(data( filter(i)+1 : filter(i+1)-1, 1));              % Relative humidity [%]
        
        
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