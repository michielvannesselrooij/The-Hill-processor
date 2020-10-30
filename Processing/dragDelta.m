function [F0, X0, dF, dFp, Xq] = dragDelta(F, X, range, ref, target)
% ------------------------------------------------------------------------
% Calculate delta drag with respect to the reference index/indeces
% By default assumes measurement scheme Ref Ref Target Ref Target Ref etc... 
%
% INPUT:
% F      : averaged force data from files
% X      : x-values for plotting, e.g. V or Re
% range  : Measurement points (velocities) to use. Default: skip last
% ref    : (optional) specify sets of references as {[a c], [d f]}
% target : (optional) specify sets of references as {[b], [e]}
% 
% OUTPUT:
% F0     : zero tared forces (cell of size F)
% dF     : absolute delta F (cell of size F)
% dFp    : relative delta F (cell of size F)
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    % Check inputs
    if sum(size(F) == size(X)) ~= 2
        if sum(size(X)) ~= 0 % allow not specifying it
            error(['The variables "F" (' num2str(size(F)) ') and "X" (' ...
                num2str(size(X)) ') should be the same size']);
        end
    end
    
    if exist('ref','var') && ~exist('target','var') || ~exist('ref','var') && exist('target','var')
        error('Specify both "target" and "ref", or neither.');
    end
    
    if exist('ref','var')
        if length(ref) ~= length(target)
            error('Variables "ref" and "target" should have the same length');
        end
    end
    
    % Default range: skip last 'null' measurement
    if ~exist('range','var')
        range = 1:length(F{1})-1;
    end
    
    % Make default X-positions if not specified
    if ~exist('X','var')
        defaultX = 1;
    elseif isempty(X)
        defaultX = 1;
    else
        defaultX = 0;
    end
    
    if defaultX
        for i = 1:length(F)
            X{i} = (1:1:length(F{i}))';
        end
    end
    
    % Make default comparison sets if not specified
    if ~exist('ref','var')    
        
        if length(F)/2 ~= floor(length(F)/2)
            t = 2:2:length(F)-1;
        else
            error('No ref/target pairs specified. Incorrect warmUp setting.');
        end
        
        for i = 1:length(t)
            ref{i}    = [t(i)-1, t(i)+1];
            target{i} =  t(i);
        end
        
    end
    
    % Trim data to range
    if exist('range','var')
        for i=1:length(F)
            F{i} = F{i}(range);
            X{i} = X{i}(range);
        end
    end
    X0 = X; % Return trimmed data
    
    % Zero all forces on first null measurement
    for i = 1:length(F)
        F{i}(isnan(F{i})) = 0;
        F0{i} = F{i}-F{i}(1);
    end
    
    for i = 1:length(target)
        
        % Define sample coordinates (target plate data)
        Xq{i} = X{target{i}};

        % Interpolate reference 'sandwhich' set on target
        n = length(ref{i});
        F_ref = zeros(length(F{1}),1);
        for j = 1:n
            F_ref = F_ref + interp1(X{ref{i}(j)}, F0{ref{i}(j)}, ...
                                        Xq{i}, 'linear', 'extrap');
        end
        
        % Average reference set
        F_ref = F_ref/n;
        
        % Calculate delta in percentage
        dF{i}  = F0{target{i}} - F_ref;
        dFp{i} = [0; dF{i}(2:end)./F_ref(2:end)*100];
        
    end

end