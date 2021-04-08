function defaultPlots(name, Cd0, Re0, dCd, dCdp, Re_target, RMSE, RMSE_X,...
    F, F_rms, F_power, p, P, X, T, Troom, pa, hum, Re, V, rho, nu, corr,...
    colors, lines, markers)
% ------------------------------------------------------------------------
% Create default output plots and save to file
%
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

%% Settings

% Ignore hotwire-only measurements
idx = find(~cellfun('isempty',F));
N   = length(idx);
if N>0

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

set(0,'DefaultLegendAutoUpdate','off');

for i = idx
    for j=1:length(dCdp{i})
        xmax(i,j) = max(Re_target{i}{j});
    end
end
xMax = max(max(xmax));

%% Summary

% Make folder for output
mkdir('results');

% Determine Reynolds number for comparison
Remax = Inf;
for i = idx
    for j=1:length(Re0{i})
        if max(Re0{i}{j}) < Remax
            Remax = max(Re0{i}{j});
        end
    end
end

% Interpolate absolute data for this Reynolds number
summary_F0 = zeros(N,1);    % Absolute force (reference)
summary_F1 = summary_F0;    % Absolute force (target)
summary_Fp0 = summary_F0;   % Pressure drag (reference)
summary_Fp1 = summary_F0;   % Pressure drag (target)
summary_Fc0 = summary_F0;   % Corrected absolute force (reference)
summary_Fc1 = summary_F0;   % Corrected absolute force (target)

i=1;
for k = idx
    for j=1:length(Re0{i})
        temp_Fc(j) = interp1(Re0{i}{j}, F{i}{j}(1:end-1)-F{i}{j}(1), Remax);
        if corr{i}{j}{2} == 0
            temp_Fp(j) = 0;
        else
            temp_Fp(j) = interp1(Re0{i}{j}, corr{i}{j}{2}(1:end-1), Remax);
        end
        temp_F(j)  = temp_Fc(j) + temp_Fp(j);
    end
    summary_F0(i)  = mean(temp_F(1:2:end));
    summary_F1(i)  = mean(temp_F(2:2:end-1));
    summary_Fp0(i) = mean(temp_Fp(1:2:end));
    summary_Fp1(i) = mean(temp_Fp(2:2:end-1));
    summary_Fc0(i) = mean(temp_Fc(1:2:end));
    summary_Fc1(i) = mean(temp_Fc(2:2:end-1));
    
    i=i+1;
end

% Create table
varNames = {'D_raw_ref', 'D_raw', 'D_p_ref', 'D_p',...
    'D_corr_ref', 'D_corr'};
T = table(summary_F0, summary_F1, summary_Fp0, summary_Fp1, summary_Fc0,...
    summary_Fc1,'VariableNames',varNames,'RowNames',name(idx));
disp(T);
writetable(T,['results' filesep 'summary.xls'],'WriteRowNames',true);

% Create charts
Dmax = max([summary_Fc0; summary_Fc1; summary_F0; summary_F1]);
i=1;
for k=idx
    
    figure;
    b=bar([summary_F0(i)-summary_Fp0(i),...
                        -summary_Fp0(i)*(summary_Fp0(i)<0),...
                         summary_Fp0(i)*(summary_Fp0(i)>0);...
           summary_F1(i)-summary_Fp1(i),...
                        -summary_Fp1(i)*(summary_Fp1(i)<0),...
                         summary_Fp1(i)*(summary_Fp1(i)>0)],...
        'stacked','FaceColor','flat');
    b(1).CData = [.1 .1 .5];
    b(2).CData = [.3 .3 .1];
    b(3).CData = [.8 .8 .8];
    ylim([0 Dmax*1.25]);
    ylabel('F [N]');
    grid minor;
    set(gca, 'XTickLabel', {'Ref', name{i}})
    legend('Corrected drag', 'Pressure correction (+)', 'Pressure correction (-)')
    
    % Save
    set(gcf,'WindowStyle','normal');
    set(gcf, 'Units', 'pixels', 'Position', [10 10 800 600]);
    saveas(gcf, ['results' filesep num2str(i) '_' name{i} '.png']);
    set(gcf,'WindowStyle','docked');
    close;
    
    i=i+1;
    
end

%% Pressure mappings

% Pressure gradient details
figure; hold on; box on; grid minor;
title('Pressure gradient along test plate');
xlabel('x [m]');
ylabel('p [Pa]');

i=1;
for k=idx
    for j=1:length(X{i})
        
        % Find pressure gradient data
        idx_grad = find(abs(X{i}{j}(:,1)) ~= 445 & X{i}{j}(:,2) == 0);
        [~, sortGrad] = sort(X{i}{j}(idx_grad,1));
        idx_grad = idx_grad(sortGrad);
         
        if length(idx_grad) >= 2
%             yPlot = P{i}{j}(:,idx_grad) ./ repmat(0.5.*rho{i}{j}.*V{i}{j}.^2, 1, length(idx_grad));
            yPlot = P{i}{j}(:,idx_grad);
            xPlot = X{i}{j}(idx_grad,1)/1000;

            for l=1:size(yPlot,1)
                
                % Plot pressure gradient data points
                if j==1 && l==2
                    h(i) = plot(xPlot, yPlot(l,:)', 'Color', c(i),...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineStyle', 'none');
                else
                    plot(xPlot, yPlot(l,:)', 'Color', c(i),...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineStyle', 'none');
                end
                
                % Fit linear pressure gradient
                pGradFit{i}{j}{l} = polyfit(xPlot, yPlot(l,:)', 1);
                pGrad{i}{j}(l)    = pGradFit{i}{j}{l}(1); % For easy plotting
                fitX              = [min(xPlot) max(xPlot)];
                fitPlot           = polyval(pGradFit{i}{j}{l}, fitX);
                
                % Plot fitted pressure gradient
                plot(fitX, fitPlot, '-', 'Color', c(i));
                
            end
        end
              
    end
    i=i+1;
end
legend(h,name(idx),'Location','South');

% Pressure gradient summary
figure; hold on; box on; grid minor;
title('Pressure gradient along test plate');
xlabel('Re [-]');
ylabel('\Deltap/\Deltax [Pa/m]');

i=1;
for k=idx
    for j=1:length(pGrad{i})

        if j==2
           h(i) = plot(Re0{i}{j}, pGrad{i}{j}(1:end-1),'Color', c(i),...
                    'Marker', m(i), 'MarkerFaceColor', c(i),'LineWidth',1.5); 
        elseif j/2 == floor(j/2) % Ref plates
            plot(Re0{i}{j}, pGrad{i}{j}(1:end-1), ':','Color', c(i),...
                    'Marker', m(i)); 
        else
            plot(Re0{i}{j}, pGrad{i}{j}(1:end-1),'Color', c(i),...
                    'Marker', m(i), 'MarkerFaceColor', c(i),'LineWidth',1.5); 
        end
        
    end
    i=i+1;
end
legend(h,name(idx),'Location','NorthWest');

% Interior pressure map
outlines{1} = [-.445, .1875;... % connector
               -.445, -.1875;...
               .445, -.1875;...
               .445, .1875;...
               -.445, .1875];
outlines{2} = [.3615, .0245;... % air bearing pockets
               .3615, .0795;...
               .4385, .0795;...
               .4385, .0245;...
               .3615, .0245];
outlines{3} = [ outlines{2}(:,1), -outlines{2}(:,2)];
outlines{4} = [-outlines{2}(:,1),  outlines{2}(:,2)];
outlines{5} = [-outlines{2}(:,1), -outlines{2}(:,2)];
outlines{6} = [.063, .010;... % sensor pin pocket
               .063, -.010;...
               .107, -.010;...
               .107, .010;...
               .063, .010];
outlines{7} = [.055, .005;... % sensor pin
               .055, -.005;...
               .065, -.005;...
               .065, .005;...
               .055, .005];
outlines{8} = [-.103, .011;... % air piston pockets
               -.103, -.011;...
               -.081, -.011;...
               -.081, .011;...
               -.103, .011];
alpha       = 30/360*2*pi;
outlines{9} = [.046 + (-0.011) *cos(alpha) -  0.011   *sin(alpha), .080 + (-0.011) *sin(alpha) +   0.011  *cos(alpha);...
               .046 + (-0.011) *cos(alpha) - (-0.011) *sin(alpha), .080 + (-0.011) *sin(alpha) + (-0.011) *cos(alpha);...
               .046 +  0.011   *cos(alpha) - (-0.011) *sin(alpha), .080 +   0.011  *sin(alpha) + (-0.011) *cos(alpha);...
               .046 +  0.011   *cos(alpha) -  0.011   *sin(alpha), .080 +   0.011  *sin(alpha) +   0.011  *cos(alpha);...
               .046 + (-0.011) *cos(alpha) -  0.011   *sin(alpha), .080 + (-0.011) *sin(alpha) +   0.011  *cos(alpha)];
outlines{10} = [outlines{9}(:,1), -outlines{9}(:,2)];

i=1;
for k=idx    
    figure;
    
    for j=1:length(X{i})
        
        for l=size(P{i}{j},1)-1
        
            subplot(length(X{i}),1,j);
            hold on; box on;
            
            if j==1
                title(['Pressure map @ V_{max} (' name(i) ')']);
            end

            % Find streamwise gap pressures
            idx_str  = find(abs(X{i}{j}(:,1)) ~= 445 & X{i}{j}(:,2) == -5);

            % Find interior pressure data
            idx_int = find(X{i}{j}(:,2) < -10);

            % Plot pressure contour
            levels = -100:1:100;
            cmap_R1 = linspace(0.2, 0.7, 128);
            cmap_R2 = linspace(0.7, 0.6, 128);
            cmap_G1 = linspace(0.2, 0.7, 128);
            cmap_G2 = linspace(0.7, 0.0, 128);
            cmap_B1 = linspace(0.6, 0.7, 128);
            cmap_B2 = linspace(0.7, 0.2, 128);
            cmap = [cmap_R1', cmap_G1', cmap_B1'; cmap_R2', cmap_G2', cmap_B2'];

            if length(idx_int) >= 3 % need at least 3 points for 2D map
                xq  = [-.445,  .445];
                zq  = [-.1875, .1875];
                xp  = X{i}{j}(idx_int,1);
                zp  = X{i}{j}(idx_int,3);
                [xq, zq] = meshgrid(xq, zq);
                vq       = griddata(xp, zp, P{i}{j}(l,idx_int), xq, zq, 'natural');
                [C, h]   = contourf(xq, zq, vq, levels);
                caxis([levels(1) levels(end)]);
                colormap(cmap);
                clabel(C,h,levels);
            end

            % Draw connector outlines
            for ii=1:length(outlines)
                plot(outlines{ii}(:,1), outlines{ii}(:,2), 'k-', 'LineWidth', 1.5);
            end

            axis equal
            xlabel('x [m]')
            ylabel('z [m]')
            xlim([-.5 .5])
            ylim([-.25 .25])
            
        end
        
    end
    i=i+1;
end          


%% Force power spectrum
i=1;
for k=idx
    
    figure; 
    
    for k=1:length(F_power{i})
        
        subplot(length(F_power{i}),1,k)
        hold on;
        box on;
        grid minor;
        ylabel('|P(f)|');
        xlim([0 50]);
        ylim([0 2000]);
        
        if k==length(F_power{i})
            xlabel('f [Hz]');
        end
        
        if k==1
            title(['Force power spectrum (' name{i} ')']);
        end

        for j=1:length(F_power{i}{k})
                plot(F_power{i}{k}{j}{1}(:,1), F_power{i}{k}{j}{1}(:,2));
                plot(F_power{i}{k}{j}{2}(:,1), F_power{i}{k}{j}{2}(:,2),...
                    'k.');
        end
    end
    
    i=i+1;
    
end

% Force power peaks vs velocity
pMax = 0;
vMax = 0;
i=1;
for l=idx
    for j=1:length(F_power{i})
        for k=1:length(F_power{i}{j})
            thres = 0.5; % axis cutoff below points with low prominance
            ptemp = F_power{i}{j}{k}{2}(F_power{i}{j}{k}{2}(:,2)>thres,1);
            pMax = max([pMax, max(ptemp)]);
            vMax = max([vMax, max(V{i}{j}(k))]);
        end
    end
    
    i=i+1;
end

%% Harmonics
i=1;
for k=idx
    
    figure;
    title(['Force signal power peaks: ' name{i} ' (valid for M-tunnel only!)']);
    hold on; box on;
%     xlabel('V [m/s]');
    xlabel('Windtunnel fan speed [rpm]');
    ylabel('f [Hz]');
%     xlim([0 vMax*1.25]);
    xlim([0 3000]);
    ylim([0 pMax*1.25]);
    
        
    % M-tunnel
    xtemp = [0, linspace(0.4,1,15)*2900, 0];
    
    % Plot rpm harmonics
    plot(xtemp(1:end-1), 0.75*xtemp(1:end-1)/60,'k:');
    plot(xtemp(1:end-1), 1.0*xtemp(1:end-1)/60,'k:');
    plot(xtemp(1:end-1), 2.0*xtemp(1:end-1)/60,'k:');
    plot(xtemp(1:end-1), 2.5*xtemp(1:end-1)/60,'k:');
    plot(xtemp(1:end-1), 5.0*xtemp(1:end-1)/60,'k:');
    plot(xtemp(1:end-1), 7.5*xtemp(1:end-1)/60,'k:');
    
    for j=1:length(F_power{i})
        
        for k=1:length(F_power{i}{j})
            power      = F_power{i}{j}{k}{2};
            pks        = power(:,1);
            prominance = power(:,2);
            [prominance, idx0] = sort(prominance);
            pks = pks(idx0);
                        
            if ~isempty(pks)
                for l=1:length(pks)
                    if prominance(l) > 0.5
                        intensity = min([1, (prominance(l)/100)^0.5]);
                        base = 0.95; % maximum fading of the point
                        dotColor = [ 1, (1-intensity)*base, (1-intensity)*base ];
                        plot(xtemp(k), pks(l), '.', 'Color', dotColor, 'MarkerSize', 10);
                        % plot(V{i}{j}(k), pks(l), '.', 'Color', dotColor, 'MarkerSize', 10);
                        % plot(xtemp(k), pks(l)/(xtemp(k)/60), '.', 'Color', dotColor, 'MarkerSize', 10);
                        
                    end
                end
            end
        end
        
    end
    
    i=i+1;
end

%% Null force shift correction

figure;
hold on;
box on;
grid minor;
title('Null force shift correction');
xlabel('Re_1 [-]')
ylabel('\Delta C_D [%]')

i=1;
for k=idx
    for j=1:length(corr{i})
        
        dF_shift_p = corr{i}{j}{1}(2:end-1) ./ (F{i}{j}(2:end-1)-F{i}{j}(1)) * 100;
        
        if j==2
            h(i) = plot(Re0{i}{j}(2:end), dF_shift_p, lines{i}, 'Color', c(i), ...
                'Marker', m(i), 'MarkerFaceColor', c(i));
        elseif j/2 ~= floor(j/2)
            plot(Re0{i}{j}(2:end), dF_shift_p, lines{i}, 'Color', c(i), ...
                'Marker', m(i));
        else
            plot(Re0{i}{j}(2:end), dF_shift_p, lines{i}, 'Color', c(i), ...
                'Marker', m(i), 'MarkerFaceColor', c(i));
        end
        
    end
    i=i+1;
end

legend(h,name(idx),'Location','EastOutside');
plot([-0.1*xMax 1.1*xMax],[0 0],'k-','LineWidth',3);
xlim([-0.1*xMax 1.1*xMax]);
ylim([ -10 10]);

% Store for dashboard
ax_null = gca;

%% Gap pressures

figure;
hold on;
grid minor;
title('Average gap pressures');
xlabel('Re_1 [-]')
ylabel('Pressure [Pa]')

i=1;
for k=idx
    for j=1:length(corr{i})
        
        if length(corr{i}{j}{3})>1   % Skip if no pressure correction
            pressureCorrectionsAvailable = 1;
            if j==2
                h1(i) = plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,2), lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineWidth', 1.5);
                h2(i) = plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,1), ':', 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineWidth', 1.5);

            elseif j==3
                h3(i) = plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,2), lines{i}, 'Color', c(i), ...
                    'Marker', m(i));
                h4(i) = plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,1), ':', 'Color', c(i), ...
                    'Marker', m(i));

            elseif j/2 == floor(j/2)
                plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,2), lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineWidth', 1.5);
                plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,1), ':', 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i), 'LineWidth', 1.5);

            else
                plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,2), lines{i}, 'Color', c(i), ...
                    'Marker', m(i));
                plot(Re0{i}{j}, corr{i}{j}{3}(1:end-1,1), ':', 'Color', c(i), ...
                    'Marker', m(i));

            end
        end
    end
    
    i=i+1;
end

% Make legend (if lines were plotted)
if exist('pressureCorrectionsAvailable','var')
    i=1;
    for k=idx
        tempName{4*(i-1)+1} = [name{i} ' - LE'];
        tempName{4*(i-1)+2} = [name{i} ' - TE'];
        tempName{4*(i-1)+3} = 'Ref - LE'; 
        tempName{4*(i-1)+4} = 'Ref - TE'; 
        h5(4*(i-1)+1) = h1(i);
        h5(4*(i-1)+2) = h2(i);
        h5(4*(i-1)+3) = h3(i);
        h5(4*(i-1)+4) = h4(i);

        i=i+1;
    end
    legend(h5,tempName,'Location','EastOutside');
end

%% Pressure drag
% (if available)

if exist('pressureCorrectionsAvailable','var')
    
    figure;
    hold on;
    box on;
    grid minor;
    title('Pressure drag');
    xlabel('Re_1 [-]')
    ylabel('F [N]')

    i=1;
    for k=idx
        for j=1:length(corr{i})

            if j==2
                h(i) = plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1), lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            elseif j/2 == floor(j/2)
                plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1), lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            else
                plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1), lines{i}, 'Color', c(i), ...
                    'Marker', m(i));
            end

        end

        i=i+1;
    end

    legend(h,name(idx),'Location','EastOutside');
    
end

%% Pressure drag contribution
% (if available)

if exist('pressureCorrectionsAvailable','var')

    figure;
    hold on;
    box on;
    grid minor;
    title('Pressure drag contribution');
    xlabel('Re_1 [-]')
    ylabel('F [%]')

    i=1;
    for k=idx
        for j=1:length(corr{i})

            Ftemp = F{i}{j}(1:end-1)-F{i}{j}(1);

            if j==2
                h(i) = plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1)./Ftemp*100, lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            elseif j/2 == floor(j/2)
                plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1)./Ftemp*100, lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            else
                plot(Re0{i}{j}, corr{i}{j}{2}(1:end-1)./Ftemp*100, lines{i}, 'Color', c(i), ...
                    'Marker', m(i));
            end

        end

        i=i+1;
    end

    legend(h,name(idx),'Location','EastOutside');

end

%% Pressure correction effect
% (if available)

if exist('pressureCorrectionsAvailable','var')

    figure;
    hold on;
    box on;
    grid minor;
    title('Pressure drag correction');
    xlabel('Re_1 [-]')
    ylabel('\Delta C_D [%]')

    % Compute pressure drag delta's
    i=1;
    for k=idx
        F_pressure = cell(size(corr{i}));
        for j=1:length(corr{i})
            F_pressure{j} = corr{i}{j}{2}(1:end-1);
        end
        [~, ~, dF_pressure{i}] = dragDelta(F_pressure, Re0{i}, 1:1:length(Re0{i}{1}));
        i=i+1;
    end

    % Calculate percentages and plot
    i=1;
    for k=idx
        for j=1:length(dF_pressure{i})

            dF_pressure_p{i}{j} = -dF_pressure{i}{j} ./ (F{i}{j}(1:end-1)-F{i}{j}(1)) * 100;

            if j==1
                h(i) = plot(Re0{i}{j}, dF_pressure_p{i}{j}, lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            else
                plot(Re0{i}{j}, dF_pressure_p{i}{j}, lines{i}, 'Color', c(i), ...
                    'Marker', m(i), 'MarkerFaceColor', c(i));
            end

        end

        i=i+1;
    end

    legend(h,name(idx),'Location','EastOutside');
    plot([-0.1*xMax 1.1*xMax],[0 0],'k-','LineWidth',3);
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([ -10 10]);

end

% Store for dashboard
ax_pcorr = gca;

%% Force

figure;
hold on;
box on;
grid minor;
title('Drag');
xlabel('Re_1 [-]')
ylabel('F [N]')

i=1;
for k = idx
    for j=1:length(Cd0{i})
        Ftemp = F{i}{j}(1:end-1)-F{i}{j}(1);
        if j==2
            h(i) = plot(Re0{i}{j},Ftemp,lines{i},'Color',c(i),'Marker',m(i),'MarkerFaceColor',c(i));
        elseif j/2 == floor(j/2)
            plot(Re0{i}{j},Ftemp,lines{i},'Color',c(i),'Marker',m(i),'MarkerFaceColor',c(i));
        else
            plot(Re0{i}{j},Ftemp,':','Color',c(i),'Marker',m(i));
        end
    end
    i=i+1;
end

legend(h,name(idx),'Location','EastOutside');

%% Cd

figure;
hold on;
box on;
grid minor;
title('Drag');
xlabel('Re_1 [-]')
ylabel('C_D [-]')

i=1;
for k = idx
    for j=1:length(Cd0{i})
        if j==2
            h(i) = plot(Re0{i}{j}(2:end),Cd0{i}{j}(2:end),lines{i},'Color',c(i),'Marker',m(i),'MarkerFaceColor',c(i));
        elseif j/2 == floor(j/2)
            plot(Re0{i}{j}(2:end),Cd0{i}{j}(2:end),lines{i},'Color',c(i),'Marker',m(i),'MarkerFaceColor',c(i));
        else
            plot(Re0{i}{j}(2:end),Cd0{i}{j}(2:end),':','Color',c(i),'Marker',m(i));
        end
    end
    i=i+1;
end

legend(h,name(idx),'Location','EastOutside');

%% Cd delta - without pressure correction
% (if available)

if exist('pressureCorrectionsAvailable','var')

    figure;
    hold on;
    box on;
    grid minor;
    title('Delta drag w/o pressure correction');
    xlabel('Re_1 [-]')
    ylabel('\Delta C_D [%]')

    i=1;
    for k = idx
        for j=1:length(dCdp{i})
            if j==1
                h(i) = plot(Re_target{i}{j}(2:end), dCdp{i}{j}(2:end) - dF_pressure_p{i}{j}(2:end),...
                lines{i}, 'Color', c(i), 'Marker', m(i),'MarkerFaceColor',c(i));
            else
                plot(Re_target{i}{j}(2:end), dCdp{i}{j}(2:end) - dF_pressure_p{i}{j}(2:end),...
                lines{i}, 'Color', c(i), 'Marker', m(i),'MarkerFaceColor',c(i));
            end

        end

        i=i+1;
    end

    legend(h,name(idx),'Location','SouthWest');
    plot([-0.1*xMax 1.1*xMax],[0 0],'k-','LineWidth',3);
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([-10 10]);

end

% Store for dashboard
ax_uncorr = gca;

%% Cd delta

figure;
hold on;
box on;
grid minor;
title('Delta drag');
xlabel('Re_1 [-]')
ylabel('\Delta C_D [%]')

i=1;
for k = idx
    for j=1:length(dCdp{i})
        if j==1
            h(i) = plot(Re_target{i}{j}(2:end), dCdp{i}{j}(2:end),...
            lines{i}, 'Color', c(i), 'Marker', m(i),'MarkerFaceColor',c(i));
        else
            plot(Re_target{i}{j}(2:end), dCdp{i}{j}(2:end),...
            lines{i}, 'Color', c(i), 'Marker', m(i),'MarkerFaceColor',c(i));
        end
        
    end
    
    i=i+1;
end

legend(h,name(idx),'Location','SouthWest');
plot([-0.1*xMax 1.1*xMax],[0 0],'k-','LineWidth',3);
xlim([-0.1*xMax 1.1*xMax]);
ylim([-10 10]);

% Store for dashboard
ax_main = gca;

%% Cd delta (averaged)

figure;
hold on;
box on;
grid minor;
title('Delta drag (averaged)');
xlabel('Re_1 [-]')
ylabel('\Delta C_D [%]')

i=1;
for k = idx
    
    % Plot average
    dFavg = zeros(size(dCdp{i}{1}));
    for j=1:length(dCdp{i})
        dFavg = dFavg + dCdp{i}{j};
    end
    dFavg = dFavg/length(dCdp{i});
    h(i) = plot(Re_target{i}{1}(2:end), dFavg(2:end), lines{i}, 'Color',...
        c(i), 'Marker', m(i), 'MarkerFaceColor', c(i), 'DisplayName', name{i});
    
    % Plot spread
%     plot(Re_target{i}{1}(2:end), dFavg(2:end)+RMSE{i}(2:end)', '--',...
%         'Color', c(i));
%     plot(Re_target{i}{1}(2:end), dFavg(2:end)-RMSE{i}(2:end)', '--',...
%         'Color', c(i));

    i=i+1;
end

legend(h,name(idx),'Location','EastOutside');
plot([-0.1*xMax 1.1*xMax],[0 0],'k-','LineWidth',3);
xlim([-0.1*xMax 1.1*xMax]);
ylim([-10 10]);

%% Delta bar chart
dF_atMax = deltaBarChart(Re_target(idx), dCdp(idx), name(idx));

%% RMSE (%)

figure;
hold on;
box on;
grid on;
title('RMSE');
xlabel('Re_1 [-]')
ylabel('\Delta C_D [%]')
ylim([0 2]);

i=1;
for k=idx
    avg = zeros(size(dCdp{i}{1}));
    for j=1:length(dCdp{i})
        avg = avg + dCdp{i}{j};
    end
    avg = avg/length(dCdp{i});
    for j=1:length(dCdp{i})
        plot(Re_target{i}{j}(2:end), abs(dCdp{i}{j}(2:end)-avg(2:end)),...
            '.:', 'Color', c(i));
    end
    
    h(i) = plot(RMSE_X{i}(2:end), RMSE{i}(2:end),lines{i}, 'Color', c(i), 'LineWidth',2);
    
    i=i+1;
end
legend(h,name(idx),'Location','NorthEast');

% Store for dashboard
ax_rmse = gca;

%% Dashboard figure
figure;
s1 = subplot(4,2,1);
    title('RMSE')
    box on;
    grid minor;
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([0 2]);
    ylabel('\Delta C_D [%]')
s2 = subplot(4,2,2);
    title('w/o Pressure correction')
    box on;
    grid minor;
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([-10 10]);
s3 = subplot(4,2,3);
    title('Null force correction')
    box on;
    grid minor;
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([-5 5]);
    xlabel('Re_1 [-]')
    ylabel('\Delta C_D [%]')
s4 = subplot(4,2,4);
    title('Pressure correction')
    box on;
    grid minor;
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([-10 10]);
    xlabel('Re_1 [-]')
    ylabel('\Delta C_D [%]')
s5 = subplot(4,2,[5:8]);
    title('Final result')
    box on;
    grid minor;
    xlim([-0.1*xMax 1.1*xMax]);
    ylim([-10 10]);
    xlabel('Re_1 [-]')

fig1 = get(ax_rmse,'children');
fig2 = get(ax_uncorr,'children');
fig3 = get(ax_null,'children');
fig4 = get(ax_pcorr,'children');
fig5 = get(ax_main,'children');

copyobj(fig1,s1);
copyobj(fig2,s2);
copyobj(fig3,s3);
copyobj(fig4,s4);
copyobj(fig5,s5);

leg_idx = [];
ii=1;
for i=1:length(fig5)
    if ~isempty(fig5(i).DisplayName)
        leg_idx(ii) = i;
        ii=ii+1;
    end
end

legend(s5,flipud(fig5(leg_idx)),'Location','SouthOutside')

end