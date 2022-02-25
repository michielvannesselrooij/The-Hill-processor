function [CI, CI_X] = calcCI(E,X)
% ------------------------------------------------------------------------
% Calculate uncertainty at 95% confidence, based on an Error input
%
% INPUT:
% E      : Cell containing equal sized vectors with the 'errors' or delta's
% X      : Cell containing equal sized vectors with query points (e.g. Re)
%
% MvN 2021 - Dimple Aerospace BV
% ------------------------------------------------------------------------

    % Check input
    if ~(size(E,1) == 1 || size(E,2) == 1) || ~iscell(E)
        error('"E" should be a one-dimensional cell variable');
    end
    
    l = zeros(length(E),2);
    for i=1:length(E)
        l(i,:) = size(E{i});
    end
    
    if sum(find(l(:,1)-mean(l(:,1)))) || sum(find(l(:,2)-mean(l(:,2))))
        error('"E" should contain equal length vectors');
    end
    
    if sum(find(l(:,2)-1)) > 0
        error('"E" should contain one-dimensional column vectors');
    end
    
    % If X not specified use integers
    if ~exist('X','var')
        X = cell(size(E));
        for i=1:length(E)
            X{i} = 1:length(E{i});
        end
    end
    
    for i=1:length(X)
        if length(X{i}) ~= length(E{i})
            error(['Size mismatch between the "X" and "E" data in cell' ...
                num2str(i)]);
        end
    end
    
    N = length(E);
    
    % Reshape the data
    E2 = zeros(N,length(E{1}));
    X2 = E2;
    for i = 1:N
        E2(i,:) = E{i};
    end
    E = E2;
    
    % Find the average x-coordinates
    for i = 1:N
        X2(i,:) = X{i};
    end
    X = X2;
    CI_X = mean(X,1);
    
    % Sample all errors E on the average x-coordinates
    for i = 1:N
        E2(i,:) = interp1(X(i,:), E(i,:), CI_X, 'linear', 'extrap');
    end
    E = E2;
    
    % Find average E
    df    = N-1;
    if df >= 1
        sigmas = tinv(0.975, df);
        CI = sigmas * std(E)/sqrt(N);
    else
        CI = NaN;
    end

end