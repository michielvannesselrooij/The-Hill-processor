function [x,y] = extractDataFromFigure
% Retrieve x-y data from current figure

fig = gcf;

dataObjs = findobj(fig,'-property','XData')
for i=1:length(dataObjs)
    x{i} = dataObjs(i).XData;
end

dataObjs = findobj(fig,'-property','YData')
for i=1:length(dataObjs)
    y{i} = dataObjs(i).YData;
end