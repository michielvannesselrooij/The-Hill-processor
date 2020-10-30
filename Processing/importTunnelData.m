function [tunnel_q, tunnel_T, tunnel_pa, tunnel_hum, tunnel_rpm] = ...
    importTunnelData(fileName, tunnel, q, T, pa, hum)
% ------------------------------------------------------------------------
% Extracts tunnel velocity and environmental data from .csv file
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    data = importdata(fileName);
    
    % Separate measurements via empty rows
    filter = [0; find(sum(abs(data),2) == 0)];
    
    % Initialize
    tunnel_T   = zeros(length(filter)-1,1);
    tunnel_q   = tunnel_T;
    tunnel_pa  = tunnel_T;
    tunnel_hum = tunnel_T;
    tunnel_rpm = tunnel_T;
    
    % Extract
    j=1;
    for i = 1:length(filter)-1
        
        if strcmp(tunnel,'M-tunnel')
        
        tunnel_q(j)     = mean(data( filter(i)+1 : filter(i+1)-1, 7));     % Dynamic pressure (q)
        tunnel_T(j)     = mean(data( filter(i)+1 : filter(i+1)-1, 8));     % Temperature test section [C]
        tunnel_pa(j)    = mean(data( filter(i)+1 : filter(i+1)-1, 9))*100; % Absolute pressure [Pa]
        tunnel_rpm(j)   = mean(data( filter(i)+1 : filter(i+1)-1, 10));    % Relative humidity [%]
        tunnel_hum      = [];
        
        % elseif strcmp(tunnel,'M-tunnel')
        
        end
        
        j=j+1;
        
    end

    % Check data and determine what to use
    if ~isempty(tunnel_q)
        if size(tunnel_q) ~= size(q)
            warning('The loaded tunnel data is not used, size mismatch.');
        end
    end
    if ~isempty(tunnel_T)
        if size(tunnel_T) ~= size(T)
            warning('The loaded tunnel data is not used, size mismatch.');
        end
    end
    if ~isempty(tunnel_pa)
        if size(tunnel_pa) ~= size(q)
            warning('The loaded tunnel data is not used, size mismatch.');
        end
    end
    if ~isempty(tunnel_hum)
        if size(tunnel_hum) ~= size(q)
            warning('The loaded tunnel data is not used, size mismatch.');
        end
    end
    
end