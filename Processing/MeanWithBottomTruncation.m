function Fmean = MeanWithBottomTruncation(Fraw, tolerance, timeout, showPlot)
% Find an improved estimate of the mean value of a signal that is truncated
% at a lower bound.
%
% INPUT
% Fraw          1D vector   Raw data time-series signal
%
% (Optional)
% tolerance     scalar      Minimal relative step size in the outcome of the iterative process
% timeout       scalar      Maximum time the iterative process may take
% showPlot      boolean     Option to show progress in a figure
%
% OUTPUT
% Fmean         scalar      Estimate of true mean of Fraw

%% Settings
if ~exist('tolerance', 'var')
    tolerance = 0.1/100;
end

if ~exist('timeout', 'var')
    timeout = 10;
end

if ~exist('showPlot', 'var')
    showPlot = false;
end

%% First estimate
Fmean = mean(Fraw);

%% Prepare plot
if showPlot
    figure;
    box on;
    hold on;
    plot(Fraw, 'b-');
    yline(Fmean, 'k')
    xlabel('Sample')
    ylabel('')
end

%% Estimate the truncation level
[~, edges, ~] = histcounts(Fraw, 100, 'Normalization', 'count');
lb = mean(edges(1:2)); % Center of first bin

%% Iterate mean estimate
ready = false;
timedout = false;
tic;

while ~(ready || timedout)
    
    previousEstimate = Fmean;

    % Estimated upper bound assuming signal symmetry
    offset = abs(lb - previousEstimate);
    ub = previousEstimate + offset; 
    
    % Clip signal at top
    Fclipped = Fraw;
    Fclipped(Fclipped > ub) = ub; 
    
    % Calculate new mean
    Fmean = mean(Fclipped); 

    % Show current estimate in plot
    if showPlot
        yline(Fmean, 'k');
    end
    
    % Stop if step is small
    stepSize = abs(1 - Fmean / previousEstimate);
    if stepSize < tolerance
        ready = true;
    end

    % Stop if time runs out
    if toc > timeout
        timedout = true;
    end
end


if showPlot
    % Show final mean value
    yline(Fmean, 'k', 'LineWidth', 2);

    % Update axis limits
    range = abs(max(Fraw) - min(Fraw));
    yMax = Fmean + 0.7*range;
    yMin = Fmean - 0.7*range;
    ylim([yMin yMax])

end