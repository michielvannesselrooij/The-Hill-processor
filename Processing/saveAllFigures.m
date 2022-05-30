function saveAllFigures(currentSize, windowSize)

%% Check inputs

% Default

if exist('currentSize','var')
    if ~islogical(currentSize) && currentSize ~= 0 && currentSize ~= 1
        error('First argument should be a boolean, if any');
    end
else
    % Do not use current size
    currentSize = false; 
end

if exist('windowSize','var')
    if numel(windowSize) ~= 2
        error('The window size should be a vector of size 1x2');
    end
else
    % Use default size
    windowSize = [800 600];
end

%% Replace / make a figure output folder
if exist('figures2','dir') == 7 
    rmdir('figures', 's');
end
mkdir figures;

%% Save figures
% Define the list of figures
figList = handle( sort( double(findall(0, 'type', 'figure') ) ));

% Save Matlab figures as .fig file
savefig(figList, ['figures' filesep 'figures'], 'compact');

% Ssave pixel images one by one
for i=1:length(figList)
    
    % Resize?
    if ~currentSize
        set(figList(i),'WindowStyle','normal');
        set(figList(i), 'Units', 'pixels', 'Position', [10; 10; windowSize(:)]');
    end
    
    saveas(figList(i), ['figures' filesep 'figure_' num2str(i) '.png']);
    
    % Dock figure again
    if ~currentSize
        set(figList(i),'WindowStyle','docked');
    end
end