function Yq = deltaBarChart(X, Y, name, Xq)
% Creates a bar chart of drag deltas at the lowest max velocity
% ------------------
% X:    cell, 2 layers
% Y:    cell of size X, 2 layers
% name: cell of size X with legend names
% Xq:   (optional) value of X for which the bar chart should be created
% ------------------

% Checks
N = length(X);
if length(Y) ~= N
    error('Provided variables should have the same length');
end


% Find average drag
for i = 1:N
    
    Yavg{i} = zeros(size(Y{i}{1}));
    for j=1:length(Y{i})
        Yavg{i} = Yavg{i} + Y{i}{j};
    end
    Yavg{i} = Yavg{i}/length(Y{i});

end


% Determine lowest max value of X
if ~exist('Xq','var')
    Xq = Inf;
    for i=N
        for j=1:length(X{i})
            if max(X{i}{j}) < Xq
                Xq = max(X{i}{j});
            end
        end
    end
end


% Interpolate at determined value
Yq = zeros(N,1);
for i=1:N
    Yq(i) = interp1(X{i}{1}, Yavg{i}, Xq, 'linear', 'extrap');
end


% Make bar chart
figure;
c = categorical(name);
b = bar(c, Yq');
ylabel('\DeltaC_D [%]');
title(['x = ' num2str(Xq, 2)]);
