function [q, dF_gap, P, X, p_gap] = pressureTranslator(p, pConfig)
% ------------------------------------------------------------------------
% Interprets pressure data from output files based on specified pressure
% channel configuration. Configuration is read from excel file
% 
% INPUT:
% p       : averaged pressures from measurement file
% pConfig : pressure channel configuration id from excel file
% 
% OUTPUT:
% q       : channel from p corresponding with the pitot delta
% dF      : calculated pressure drag correction
% P       : mapped pressures
% X       : position vectors of P
% p_gap   : average pressure in gap
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------

% Check if pressure data is given
if p==0
    q      = 0;
    dF_gap = 0;
    P      = 0;
    X      = 0;
    p_gap  = 0;
else

% Find mastersheet location
current    = cd;                                             % Current path
baseFolder = [filesep 'Analyses'];                           % Analysis folder
idx        = strfind(current,baseFolder)-length(baseFolder); % Trim path idx
base       = current(1:idx+length(baseFolder)-1);            % Trim path

file   = [base filesep 'Measurements' filesep 'mastersheet.xlsx'];

% Load pressure configurations from Mastersheet
[~, ~, layouts] = xlsread(file, 'Pressure layouts');
config = 1;
for i=1:size(layouts,2)
    if strcmp(layouts{2,i}, pConfig)
        config = cell2mat(layouts(5:end,i));
        break;
    end
end

% Zero the pressures
p = p-p(1,:);

% Correct for residual pressures at end of measurement
p_shift  = p(end,:);
p(end,:) = zeros(size(p(end,:)));
p        = p - repmat(p_shift,size(p,1),1).*abs(p)./max(abs(p));

% Extract pitot pressure
if isnan(config(1))
    q = [];
else
    q = p(:,config(1));
end

% Map pressures
channels = config(~isnan(config));
idx      = find(~isnan(config));    % Find active pressure taps

if length(idx)<2
    dF_gap = 0;
    p_gap  = 0;
    P      = [];
    X      = [];
    
else
    for i=2:length(idx)
        X(i-1,1) = layouts{4+idx(i),5}; % X location
        X(i-1,2) = layouts{4+idx(i),6}; % Y location
        X(i-1,3) = layouts{4+idx(i),7}; % Z location
        P(:,i-1) = p(:,channels(i));    % Pressure @ all V
    end

    % Extract original LE/TE pressure tap locations for plotting
    yLE_plot = X(X(:,1)==-445,2)/1000;
    zLE_plot = X(X(:,1)==-445,3)/1000;
    yTE_plot = X(X(:,1)==445,2)/1000;
    zTE_plot = X(X(:,1)==445,3)/1000;

    % --------------------------------------------------------------------
    % Calculate pressure correction --------------------------------------
    % --------------------------------------------------------------------

    % Prepare
    idx = idx(2:end)-1;                    % Remove pitot from channels
    yq  = linspace(-10,0,5)/1000;          % Interpolation points y
    zq  = linspace(-187.5,187.5,100)/1000; % Interpolation points z

    % Vertical pressure gradient (check if the TE center taps are used)
    % --------------------------------------------------------------------
    p_vert = cell(size(p,1),1);
    lower = find(idx==19);
    upper = find(idx==20);

    if length(lower) + length(upper) == 2

        p1    = P(:,lower);     % lower pressure
        p2    = P(:,upper);     % upper pressure    
        pmean = mean([p1 p2],2);

        % Calculate vertical pressure distribution for each velocity
        for i = 1:size(p,1) 
            p_vert{i} = interp1([-7.5, -2.5]/1000,...
                [p1(i)./pmean(i), p2(i)./pmean(i)], yq, 'linear', 'extrap');
        end

        % Replace center measurements with their average
        P(:,[lower upper]) = [];
        P = [P pmean];
        X(upper,2) = -5;
        X(lower,:) = [];

    % Vertical pressure distribution based on one sample only
    elseif length(lower) + length(upper) == 1

        warning(['Only one of the center TE pressure taps are used!'...
            'Results may be invalid.']);

        % Interpolate center location
        idx_temp           = find(X(:,2)==-5 & X(:,1)==445);
        z_temp             = X(idx_temp,3);
        p_temp             = P(:,idx_temp);
        [z_temp, idx_sort] = sort(z_temp);
        p_temp             = p_temp(:,idx_sort);

        for i=1:size(p,1)
            p_estimate(i) = interp1(z_temp, p_temp(i,:), 0);
        end

        if length(lower) == 1
            p1 = p(:,lower);  % lower pressure
            p2 = p_estimate;
            y1 = -7.5;
            y2 = -5;

            % Replace with center estimate for further processing
            P(:,lower) = p_estimate;
            X(lower,2) = -5;

        else
            p1 = p_estimate;
            p2 = p(:,upper);  % upper pressure
            y1 = -5;
            y2 = -2.5;

            % Replace with center estimate for further processing
            P(:,upper) = p_estimate;
            X(upper,2) = -5;

        end

        % Calculate vertical pressure distribution for each velocity
        pmean  = mean([p1 p2],2);
        for i = 1:size(p,1) 
            p_vert{i} = interp1([y1, y2]/1000,...
                [p1(i)./pmean(i), p2(i)./pmean(i)], yq, 'linear', 'extrap');
        end


    else
        % No vertical gradient adjustment
        for i=1:length(p_vert)
            p_vert{i} = ones(1,length(yq));
        end

    end

    % Horizontal pressure gradient
    % --------------------------------------------------------------------
    LE  = find(X(:,2)==-5 & X(:,1)==-445);        % Select LE tap indices
    TE  = find(X(:,2)==-5 & X(:,1)==445);         % Select TE tap indices

    % Leading edge
    [z_LE, idx_temp] = sort(X(LE,3)');
    p_LE = P(:,LE(idx_temp));

    for i=1:size(p_LE,1)

        p_int_LE{i} = interp1(z_LE/1000, p_LE(i,:), zq, 'spline', 'extrap');

    %     p_int_LE{i} = repmat(p_int_LE{i}, length(yq), 1)...
    %                     .* repmat(p_vert{i}', 1, length(zq));

        % Don't use vertical gradient for LE
        p_int_LE{i} = repmat(p_int_LE{i}, length(yq), 1); 

    end

    % Trailing edge
    [z_TE, idx_temp] = sort(X(TE,3)');
    p_TE = P(:,TE(idx_temp));

    for i=1:size(p_TE,1)

        p_int_TE{i} = interp1(z_TE/1000, p_TE(i,:), zq, 'spline', 'extrap');

    %     p_int_TE{i} = repmat(p_int_TE{i}, length(yq), 1)...
    %                     .* repmat(p_vert{i}', 1, length(zq));
        p_int_TE{i} = repmat(p_int_TE{i}, length(yq), 1)...
                        + repmat(p_vert{i}', 1, length(zq));                

    end

    % Pressure delta & drag correction
    dF_gap = zeros(length(p_int_TE),1);
    dp     = cell(length(p_int_TE),1);
    for i=1:length(p_int_TE)

        dp{i}     = p_int_LE{i} - p_int_TE{i};

        Aq        = ( (zq(2)-zq(1)) * (yq(2)-yq(1)) );
        dF_gap(i) = sum(sum(Aq.*dp{i}));

    end

    dF_gap(isnan(dF_gap)) = 0;

    % Plot gap pressures
    figure; 
    levels = -100:10:100;

    cmap_R1 = linspace(0.2, 0.3, 128);
    cmap_G1 = linspace(0.2, 0.3, 128);
    cmap_B1 = linspace(0.6, 0.3, 128);
    cmap_R2 = linspace(0.3, 0.6, 128);
    cmap_G2 = linspace(0.3, 0.0, 128);
    cmap_B2 = linspace(0.3, 0.2, 128);
    cmap = [cmap_R1', cmap_G1', cmap_B1'; cmap_R2', cmap_G2', cmap_B2'];

    subplot(3,1,1); box on; ylabel('y [m]'); hold on;
    [C1, h1] = contourf(zq, yq', p_int_LE{end-1}, levels);
    scatter(zLE_plot, yLE_plot, 10, 'w');
    title('LE gap p [Pa]'); 
    caxis([levels(1) levels(end)]);
    colormap(cmap);

    subplot(3,1,2); box on; ylabel('y [m]'); hold on;
    [C2, h2] = contourf(zq, yq, p_int_TE{end-1}, levels);
    scatter(zTE_plot, yTE_plot, 10, 'w');
    title('TE gap p [Pa]'); 
    caxis([levels(1) levels(end)]);
    colormap(cmap);

    subplot(3,1,3); box on; ylabel('y [m]'); xlabel('z [m]'); hold on;
    [C3, h3] = contourf(zq, yq, dp{end-1}, levels);
    title('\Delta p [Pa]'); 
    caxis([levels(1) levels(end)]);
    colormap(cmap);

    clabel(C1,h1,levels);
    clabel(C2,h2,levels);
    clabel(C3,h3,levels);

    % --------------------------------------------------------------------
    % Substitute interpolated data in output -----------------------------
    % --------------------------------------------------------------------

    P(:,[LE; TE]) = [];
    X([LE; TE],:) = [];

    P2 = zeros(size(P,1), length(zq)*length(yq));
    P3 = P2;
    X2 = zeros(length(zq)*length(yq),3);
    X3 = X2;

    n1 = size(p_int_TE{i},1);
    n2 = size(p_int_TE{i},2);
    for j=1:n1
        for k=1:n2

            for i=1:length(p_int_TE)
                P2(i,(j-1)*n2+k) = p_int_TE{i}(j,k);
            end
            X2((j-1)*n2+k,:) = [445, yq(j), zq(k)];

        end
    end

    for j=1:n1
        for k=1:n2

            for i=1:length(p_int_LE)
                P3(i,(j-1)*n2+k) = p_int_LE{i}(j,k);
            end
            X3((j-1)*n2+k,:) = [-445, yq(j), zq(k)];

        end
    end

    % Pressures and coordinates
    P = [P, P2, P3];
    X = [X; X2; X3];

    % Average gap pressures
    p_gap = [mean(P2,2), mean(P3,2)];
end

end