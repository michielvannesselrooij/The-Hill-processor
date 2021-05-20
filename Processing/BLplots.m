function BLplots(name, y, u, u_rms, u_power, nu, ut, y0, k, B, PI, d, d_star,...
            theta, H, up_model, yp_model, V, dCdp, colors, lines, markers)

% DEFAULT BOUNDARY LAYER PLOTS
% ----------------------------------------------------------
% Accepts none, one or two hotwire data sets per 'name'
%   0 data sets: skipped in all plots
%   1 data set:  Single line, skipped in delta plots
%   2 data sets: 1st measurement is reference, 2nd is target
% ----------------------------------------------------------


%% Analyze number of measurements
N = length(name);
for i=1:N
    n(i) = length(u{i});
end

n_empty  = find(n==0); % number of entries with hotwire data
n_single = find(n==1); % number of single hotwire measurement
n_double = find(n==2); % number of hotwire delta sets
n_limit  = find(n>2);  % number of incorrect number of inputs

% Warning if more than 2 hotwire sets are given for one entry
if ~isempty(n_limit)
    warning('Too many hotwire samples specified!');
end

% Stop if no hotwire data found
if isempty(n_single) && isempty(n_double)
    return;
end

%% Preparations


%% Styling settings

% Determine colors
c = '';
for i=1:N
    if isnan(colors{i})
        c   = 'krbmgckrbmgc'; % Default
        break;
    else
        c = [c, colors{i}];
    end
end

% Determine markers
m = '';
for i=1:N
    if isnan(markers{i})
        m   = 's<d>os<d>os<'; % Default
        break;
    else
        m = [m, markers{i}];
    end
end

% Determine line styling
for i=1:N
    if isnan(lines{i})
        for j=1:N
            lines{j} = '-'; % Default
        end
        break;
    end
end

c2 = distinguishable_colors(N);    % Colors  (consolidations)

idx          = 1:N;
idx_single   = idx(n_single);
idx_double   = idx(n_double);
idx(n_empty) = [];
N            = length(idx);

%% Diagnostic plot (Alfredson, 2010)
for i=idx
    figure; hold on; box on; grid on;
    xlim([0 1]);
    ylim([0 0.12]);
    xlabel('u/U [-]');
    ylabel("u'/U [-]");
    title(['Diagnostic plot - ' name(i)]);
    
    plot([0 0.2], [0 0.2*0.4], 'k-.');         % Near-wall gradient line
    plot([0.2 1], [0.2 1]*0.157+0.047, 'k--'); % Estimated maximum line
        
    for j=1:length(u{i})
        h(j) = plot(u{i}{j}/u{i}{j}(end), u_rms{i}{j}/u{i}{j}(end),...
            'o','Color',c(j));
        
        % Estimated maximum
        Umax = 10.6*ut{i}(j)/u{i}{j}(end);
        plot(Umax, Umax*0.157+0.047,...
            'o','Color',c(j),'MarkerFaceColor',c(j))
    end
    
    if length(u{i})==2
        legend(h,'Ref.','Target','Location','South');
    end
end

%% Key boundary layer properties (absolute)
figure;

j=1;
for i=idx

    if ~isempty(find(idx_single==i)) % in case of single measurement
        plot_prop{1}(j,:)        = [0, ut{i}];
        plot_prop{2}(j,:)        = [0, d{i}*1e3];
        plot_prop{3}(j,:)        = [0, d_star{i}*1e3];
        plot_prop{4}(j,:)        = [0, theta{i}*1e3];
        plot_prop{5}(j,:)        = [0, H{i}];
        plot_prop{6}(j,:)        = [0, PI{i}];
        plot_prop{7}(j,:)        = [0, k{i}];
        plot_prop{8}(j,:)        = [0, B{i}];
        
    else
        plot_prop{1}(j,:)        = ut{i};
        plot_prop{2}(j,:)        = d{i}*1e3;
        plot_prop{3}(j,:)        = d_star{i}*1e3;
        plot_prop{4}(j,:)        = theta{i}*1e3;
        plot_prop{5}(j,:)        = H{i};
        plot_prop{6}(j,:)        = PI{i};
        plot_prop{7}(j,:)        = k{i};
        plot_prop{8}(j,:)        = B{i};
    end
    
    j=j+1;
    
end

label_x = { ...
    cell(size(idx)),...
    cell(size(idx)),...
    cell(size(idx)),...
    cell(size(idx)),...
    cell(size(idx)),...
    cell(size(idx)),...
    cell(size(idx)),...
    name,...
    };

label_y = { ...
    'u_\tau [m/s]',...
    '\delta [mm]',...
    '\delta^* [mm]',...
    '\Theta [mm]',...
    'H [-]',...
    '\Pi [-]',...
    '\kappa [-]',...
    'B [-]',...
    };

subTit = { ...
    'u_\tau',...
    '\delta',...
    '\delta^*',...
    '\Theta',...
    'H',...
    '\Pi',...
    '\kappa',...
    'B',...
    };

for i=1:length(plot_prop)
    subplot(length(plot_prop),1,i); box on; 
    ytemp = plot_prop{i};
    if isrow(ytemp)
        ytemp = vertcat(ytemp,nan(size(ytemp)));
    end
    b = bar(ytemp, 'grouped');
    b(1).FaceColor = c(1);
    b(2).FaceColor = c(2);
    if isrow(plot_prop{i})
        xlim([0.5 1.5])
    end
    set(gca, 'XTickLabel', label_x{i})
    ylabel(label_y{i})
    title(subTit{i})
    ylim([0.95*min(min(plot_prop{i})) 1.05*max(max(plot_prop{i}))]);
    for j=1:size(plot_prop{i},2)
        text(b(j).XData+b(j).XOffset, b(j).YData, ...
                    num2str(plot_prop{i}(:,j),'%0.3f'), ...
                   'HorizontalAlignment','center', ...
                   'VerticalAlignment','bottom')
    end
end

%% Property delta's

if length(idx_double)>1
    
    % Key boundary layer properties (Delta [%])
    figure; box on; grid minor;
    
    d_plot_prop = [];
    for i=1:length(plot_prop)
        d_plot_prop = [d_plot_prop, (plot_prop{i}(:,2:end) ./ ...
            repmat(plot_prop{i}(:,1),1,size(plot_prop{i},2)-1)-1)*100];
    end
    
    emptyRows = find((sum(d_plot_prop,2)==Inf));
    d_plot_prop(emptyRows,:) = [];
    
    b = bar(d_plot_prop);
    cBar = winter(size(d_plot_prop,2));
    for i=1:size(d_plot_prop,2)
        b(i).FaceColor = cBar(i,:);
    end
    set(gca, 'XTickLabel', name(idx_double))
    ylabel('[%]')
    title('Changes to boundary layer properties')
    legend(subTit,'Orientation', 'Horizontal', 'Location', 'South',...
        'FontSize',14);
    
    limitY = ceil(max(abs([ ...
        min(min(d_plot_prop)) max(max(d_plot_prop)) ]))/5)*5+5;
    ylim([ -limitY limitY ]);
    
    for i=1:size(d_plot_prop,2)
        for j=1:size(d_plot_prop,1)
            if d_plot_prop(j,i) >= 0
                text(b(i).XData(j)+b(i).XOffset, b(i).YData(j), ...
                        num2str(d_plot_prop(j,i),'%0.1f'), ...
                       'HorizontalAlignment','center', ...
                       'VerticalAlignment','bottom')
            else
                text(b(i).XData(j)+b(i).XOffset, b(i).YData(j), ...
                        num2str(d_plot_prop(j,i),'%0.1f'), ...
                       'HorizontalAlignment','center', ...
                       'VerticalAlignment','top')
            end
        end
    end

    % Boundary layer property correlation with Cd if available
    k=1;
    for i=idx_double
        
        if ~isempty(V{i})

            % Determine drag delta values at BL measurement velocities
            Vdrag  = [];
            Cddrag = [];
            
            for j=1:length(V{i})
                Vdrag = [Vdrag, V{i}{j}(1:end-1)];
            end
            Vdrag = mean(Vdrag, 2);

            for j=1:length(dCdp{i})
                Cddrag = [Cddrag, dCdp{i}{j}];
            end
            Cddrag = mean(Cddrag, 2);

            % Find drag delta at velocity of BL measurements (average)
            dragDelta(k) = 0;
            for j=1:length(u{i})
                dragDelta(k) = dragDelta(k) + ...
                    interp1(Vdrag, Cddrag, u{i}{j}(end));
            end
            dragDelta(k) = dragDelta(k)/length(u{i});
            
            % Add BL properties to plot variable
            d_plot_prop2(k,:) = d_plot_prop(k,:);
            
            k=k+1;
        end

    end

    figure;
    nPlots = size(d_plot_prop2,2);
    for i=1:nPlots
        subplot(nPlots,1,i); box on;
        plot( d_plot_prop2(:,i), dragDelta, 'kd');
        title(subTit{i})
        xlim([ -limitY limitY ])
        ylim([ min(dragDelta)-0.5 max(dragDelta)+0.5 ])
        if i==nPlots
            xlabel('Change in property [%]')
        end
        ylabel('\Delta C_D [%]')
    end

end 

%% Dimensional profile (mean)
figure;
for i=idx
    
    subplot(1,N,find(idx==i));
    hold on; grid on; box on;
    
    maxX = 0; maxY = 0;
    for j=1:length(y{i})
        if isempty(find(idx_double==i))
            h = plot(u{i}{j}, y{i}{j}*1e3, 's-','Color',c(j+1));
        else
            h(j) = plot(u{i}{j}, y{i}{j}*1e3, 's-','Color',c(j));
        end
        maxX = max([maxX, max(u{i}{j})]);
        maxY = max([maxY, max(y{i}{j})]);
    end
    
    if i==1
        ylabel('y [mm]');
    end
    
    if ~isempty(find(idx_double==i))
        legend(h, 'Ref.','Target','Location','NorthWest');
    end
    
    xlabel('u [m/s]');
    title(name{i});
    xlim([0 maxX*1.1]);
    ylim([0 maxY*1.1*1e3]);
    
end

%% Dimensional profile (turbulence)
load(['..' filesep '..' filesep 'Processing' filesep 'data_sets'...
    filesep 'klebanoff.mat'],'kleb_u','kleby_u');

figure;
for i=idx
    
    subplot(1,N,find(idx==i));
    hold on; grid on; box on;
    
    maxX = 0; maxY = 0;
    for j=1:length(y{i})
        if isempty(find(idx_double==i))
            h = plot(u_rms{i}{j}/u{i}{j}(end), y{i}{j}*1e3, 's-','Color',c(j+1));
            plot(kleb_u,kleby_u*d{i}(j)*1e3,':','Color',c(j+1)); % Klebanoff profile
        else
            h(j) = plot(u_rms{i}{j}/u{i}{j}(end), y{i}{j}*1e3, 's-','Color',c(j));
            plot(kleb_u,kleby_u*d{i}(j)*1e3,':','Color',c(j)); % Klebanoff profile
        end
        maxX = max([maxX, max(u_rms{i}{j}/max(u{i}{j}))]);
        maxY = max([maxY, max(y{i}{j})]);
    end
    
    if i==1
        ylabel('y [mm]');
    end
    
    xlabel("u'/U [-]");
    
    if ~isempty(find(idx_double==i))
        legend(h, 'Ref.','Target');
    end
    
    title(name{i});
    xlim([0 maxX*1.1]);
    ylim([0 maxY*1.1*1e3]);
    
end

%% Velocity profile (consolidated)

figure; 
hold on; 
grid on; 
box on; 
maxX = 0; 
maxY = 0; 
ii=1;
h = [];

for i=idx
    
    for j=1:length(y{i})
        
        
        if length(y{i}) == 1 || (length(y{i}) == 2 && j == 2)
            
            h(ii) = plot(u{i}{j}/u{i}{j}(end), y{i}{j}*1e3, '-',...
            'Color',c2(i,:),'LineWidth',2);
            ii=ii+1;
            
        else
            
            plot(u{i}{j}/u{i}{j}(end), y{i}{j}*1e3, '-',...
            'Color',c2(i,:),'LineWidth',2);
        
        end
        
        maxX = max([maxX, max(u{i}{j}/u{i}{j}(end))]);
        maxY = max([maxY, max(y{i}{j})]);
    end
    
end

xlabel("u/U [-]");
ylabel('y [mm]');
legend(h,name(idx),'Location','NorthWest');
title('Velocity profile (consolidated)');
xlim([0 maxX*1.1]);
ylim([0 maxY*1.1*1e3]);

%% Turbulence intensity (consolidated)

figure; 
hold on; 
grid on; 
box on; 
maxX = 0; 
maxY = 0; 
ii=1;
h = [];

for i=idx
    for j=1:length(y{i})
        
        if length(y{i}) == 1 || (length(y{i}) == 2 && j == 2)
            h(ii) = plot(u_rms{i}{j}/u{i}{j}(end), y{i}{j}*1e3, '-',...
                'Color', c2(i,:),'LineWidth',2);
            ii=ii+1;
        else
            plot(u_rms{i}{j}/u{i}{j}(end), y{i}{j}*1e3, '-',...
                'Color', c2(i,:),'LineWidth',2);
        end
        
        maxX = max([maxX, max(u_rms{i}{j}/u{i}{j}(end))]);
        maxY = max([maxY, max(y{i}{j})]);
    end
    
end

xlabel("u'/U [-]");
ylabel('y [mm]');
legend(h,name(idx),'Location','NorthEast');
title('Turbulence intensity (consolidated)');
xlim([0 maxX*1.1]);
ylim([0 maxY*1.1*1e3]);

%% Power spectra
for i=idx
    
    % Determine interpolation range for each BL
    for j=1:length(u_power{i})
        for k=1:length(u_power{i}{j})
            f_min{j}(k) = min(u_power{i}{j}{k}(:,1));
            f_max{j}(k) = max(u_power{i}{j}{k}(:,1));
        end
    end
    for j=1:length(u_power{i})
        fq_min(j) = max(f_min{j});
        fq_max(j) = min(f_max{j});
    end
    fq = linspace(max(fq_min), min(fq_max), 200);
    
    % Interpolate data
    for j=1:length(u_power{i})
        for k=1:length(u_power{i}{j})
            powerq{i}{j}(k,:) = abs(interp1(u_power{i}{j}{k}(:,1),...
                                u_power{i}{j}{k}(:,2),fq,'linear'));
        end
    end
    
    figure;
    levels = linspace(0,10,200);
    
    if isempty(find(idx_double==i))
        
        [C, h] = contourf(fq, y{i}{1}, powerq{i}{1}, levels);
        set(h,'LineColor','none')
        caxis([levels(1) levels(end)]);
        colormap('default');
        xlabel('f [Hz]'); ylabel('y [m]');
        
    else
        n = length(powerq{i});
        for j=1:n
        
            subplot(1,n,j);
            [C, h] = contourf(fq, y{i}{j}, powerq{i}{j}, levels);
            set(h,'LineColor','none')
            caxis([levels(1) levels(end)]);
            colormap('default');
            xlabel('f [Hz]'); ylabel('y [m]');
        
        end
    end
end

%% Non-dimensional profile

if isempty(idx_double)
    consolidate = true;
else
    consolidate = false;
end

if consolidate
    figure; hold on; grid on; box on;
    h=[];

    maxX = 0; maxY = 0;
    for i=1:length(y)
        for j=1
            yp{i}{j} = y{i}{j}*ut{i}(j)/nu{i}(j);
            up{i}{j} = u{i}{j}/ut{i}(j);
            h(i) = plot(yp{i}{j}, up{i}{j}, 'Color', c(i), 'Marker', '.',...
                'LineStyle', 'none');
            
            h2 = plot(yp_model{i}{j}, up_model{i}{j}, 'k-');

            maxX = max([maxX, max(yp{i}{j})]);
            maxY = max([maxY, max(up{i}{j})]);
        end
    end

    legend([h, h2], [name, 'model'], 'Location', 'NorthWest');

    set(gca,'xscale','log');
    xlim([1 maxX*3]);
    ylim([0 maxY*1.25]);
    xlabel('y^+');
    ylabel('u^+');
        
else

    figure;
    for i=idx

        subplot(N,1,find(idx==i));
        hold on; grid on; box on;

        maxX = 0; maxY = 0;
        for j=1:length(y{i})

            if isempty(find(idx_double==i))
                yp{i}{j} = y{i}{j}*ut{i}(j)/nu{i}(j);
                up{i}{j} = u{i}{j}/ut{i}(j);
                plot(yp_model{i}{j}, up_model{i}{j}, 'k-');
                h = plot(yp{i}{j}, up{i}{j}, 's', 'Color', c(j+1));
            else
                yp{i}{j} = y{i}{j}*ut{i}(j)/nu{i}(j);
                up{i}{j} = u{i}{j}/ut{i}(j);
                plot(yp_model{i}{j}, up_model{i}{j}, 'k-');
                h(j) = plot(yp{i}{j}, up{i}{j}, 's', 'Color', c(j));
            end

            maxX = max([maxX, max(yp{i}{j})]);
            maxY = max([maxY, max(up{i}{j})]);
        end

        if ~isempty(find(idx_double==i))
            legend(h, 'Ref.','Target', 'Location', 'NorthWest');
        end

        set(gca,'xscale','log');
        xlim([1 maxX*3]);
        ylim([0 maxY*1.25]);
        xlabel('y^+');
        ylabel('u^+');
        title(name{i});
    end
end