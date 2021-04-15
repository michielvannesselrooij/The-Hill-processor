function [q, dF_gap, dF_int, P, X, p_gap] = pressureTranslator(p, pConfig)
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
    dF_int = 0;
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
p_shift                 = p(end,:);
p_shift(isnan(p_shift)) = 0;
p(end,:)                = zeros(size(p(end,:)));
% p                       = p - repmat(p_shift, size(p,1), 1) .* abs(p)./max(abs(p));

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
    dF_int = 0;
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
    % Calculate gap pressure correction ----------------------------------
    % --------------------------------------------------------------------

    % Prepare
    idx = idx(2:end)-1;                    % Remove pitot from channels
    yq  = linspace(-10,0,5)/1000;          % Interpolation points y
    zq  = linspace(-187.5,187.5,100)/1000; % Interpolation points z    
    n1  = length(yq);
    n2  = length(zq);

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
            if pmean(i) == 0
                p_vert{i} = ones(1,n1);
            else
                p_vert{i} = interp1([-7.5, -2.5]/1000,...
                    [p1(i)./pmean(i), p2(i)./pmean(i)], yq, 'linear', 'extrap');
            end
        end

        % Replace center measurements with their average
        P(:,upper) = pmean;
        P(:,lower) = [];
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
        p_int_LE{i} = repmat(p_int_LE{i}, length(yq), 1)...
                        .* repmat(p_vert{i}', 1, length(zq));
    end

    % Trailing edge
    [z_TE, idx_temp] = sort(X(TE,3)');
    p_TE = P(:,TE(idx_temp));

    for i=1:size(p_TE,1)
        p_int_TE{i} = interp1(z_TE/1000, p_TE(i,:), zq, 'spline', 'extrap');
        p_int_TE{i} = repmat(p_int_TE{i}, length(yq), 1)...
                        .* repmat(p_vert{i}', 1, length(zq));
    end

    % Pressure delta & drag correction
    % --------------------------------------------------------------------
    dF_gap = zeros(length(p_int_TE),1);
    dp     = cell(length(p_int_TE),1);
    for i=1:length(p_int_TE)

        dp{i}     = p_int_LE{i} - p_int_TE{i};

        Aq        = ( (zq(2)-zq(1)) * (yq(2)-yq(1)) );
        dF_gap(i) = sum(sum(Aq.*dp{i}));

    end

    dF_gap(isnan(dF_gap)) = 0;
    
    
    % --------------------------------------------------------------------
    % Calculate interior pressure correction -----------------------------
    % --------------------------------------------------------------------
    
    % Find interior pressure data
    dx      = 1e-3; % discretization [m]
       
    % Specify vertical surfaces
    faceParams = [...
    %    x_min   x_max    z_min   z_max  h     x-normal
         .3615,  .3615,  .0245,  .0795, .004,  1;... % air bearing pockets (inside)
         .3651,  .3651, -.0245, -.0795, .004,  1;...
        -.3651, -.3651, -.0245, -.0795, .004, -1;...
        -.3651, -.3651,  .0245,  .0795, .004, -1;...
         .4385,  .4385,  .0245,  .0795, .004, -1;... % air bearing pockets (outside)
         .4385,  .4385, -.0245, -.0795, .004, -1;...
        -.4385, -.4385, -.0245, -.0795, .004,  1;...
        -.4385, -.4385,  .0245,  .0795, .004,  1;...
         .0630,  .0630, -.0100,  .0100, .004,  1;... % Sensor pin pocket
         .1070,  .1070, -.0100,  .0100, .004, -1;...
         .0550,  .0550, -.0050,  .0050, .010, -1;... % Sensor pin
         .0650,  .0650, -.0050,  .0050, .010,  1;...
        -.1030, -.1030, -.0110,  .0110, .004,  1;... % Air piston pockets (center)
        -.0810, -.0810, -.0110,  .0110, .004, -1;...
         .0310,  .0420, -.0840, -.0650, .004,  cos(30/360*2*pi);... % Air piston pockets (angled)
         .0310,  .0500, -.0840, -.0950, .004,  sin(30/360*2*pi);...
         .0610,  .0420, -.0760, -.0650, .004, -cos(30/360*2*pi);...
         .0610,  .0500, -.0760, -.0950, .004, -sin(30/360*2*pi);...
         .0310,  .0420,  .0840,  .0650, .004,  cos(30/360*2*pi);...
         .0310,  .0500,  .0840,  .0950, .004,  sin(30/360*2*pi);...
         .0610,  .0420,  .0760,  .0650, .004, -cos(30/360*2*pi);...
         .0610,  .0500,  .0760,  .0950, .004, -sin(30/360*2*pi);...
        ];

    face = {}; % [x, z, h, x-normal]
    i=1;
    for k=1:size(faceParams,1)

        % Streamwise-facing surfaces
        if abs( faceParams(k,6) ) == 1
            if faceParams(k,4) > faceParams(k,3)
                z = faceParams(k,3) : dx : faceParams(k,4);
            else
                z = fliplr( faceParams(k,4) : dx : faceParams(k,3) );
            end

            face{i}(:,1) = faceParams(k,1) * ones(size(z));
            face{i}(:,2) = z;
            face{i}(:,3) = faceParams(k,5) * ones(size(z));
            face{i}(:,4) = faceParams(k,6) * ones(size(z));
        
        % Other surfaces
        else
            n = floor( sqrt( (faceParams(k,2)-faceParams(k,1))^2 ...
                            +(faceParams(k,4)-faceParams(k,3))^2 )/dx );
            dx_x = abs(faceParams(k,2)-faceParams(k,1)) / n;
            dx_z = abs(faceParams(k,4)-faceParams(k,3)) / n;
            
            if faceParams(k,2) > faceParams(k,1)
                face{i}(:,1) = faceParams(k,1) : dx_x : faceParams(k,2);
            else
                face{i}(:,1) = fliplr( faceParams(k,2) : dx_x : faceParams(k,1) );
            end
            
            if faceParams(k,4) > faceParams(k,3)
                face{i}(:,2) = faceParams(k,3) : dx_z : faceParams(k,4);
            else
                face{i}(:,2) = fliplr( faceParams(k,4) : dx_z : faceParams(k,3) );
            end
            
            face{i}(:,3) = faceParams(k,5) * ones(n+1,1);
            face{i}(:,4) = faceParams(k,6) * ones(n+1,1);
        end
        
        % Discretized face centers
        face_center{i}(:,1) = (face{i}(1:end-1,1) + face{i}(2:end,1)) / 2;
        face_center{i}(:,2)  = (face{i}(1:end-1,2) + face{i}(2:end,2)) / 2;
        face_center{i}(:,3)  = face{i}(1:end-1,3);
        face_center{i}(:,4)  = face{i}(1:end-1,4);

        i=i+1;
    end
    
%     % Show configuration (optional)
%     figure; hold on; box on; axis equal;
%     title('Interior faces used in correction')
%     xlim([-.5 .5]);
%     ylim([-.2 ,.2]);
%     for i=1:length(faces)
%         
%         if faces{i}(1,4) > 0
%             cDir = 'r';
%             
%         else
%             cDir = 'b';
%         end
%         
%         plot(face{i}(:,1),face{i}(:,2),'-', 'Color', cDir)
%         plot(face_center{i}(:,1),face_center{i}(:,2),'k.', 'Color', cDir)
%         
%     end
    
    % Calculate pressure correction
    idx_int = find(X(:,2) < -10);
    dF_int = zeros(size(P,1),1);
    x_int   = X(idx_int,1)';
    z_int   = X(idx_int,3);
    p_int   = P(:,idx_int);
%     [x_int, z_int] = meshgrid(x_int, z_int);
    
    if length(idx_int) >= 3 % need at least 3 points for 2D map

        for i=1:length(face)

            for j=1:size(P,1)

                % Pressure on the face for each velocity
                pFace = griddata(x_int, z_int, p_int(j,:), face_center{i}(:,1), face_center{i}(:,2));

                % Add drag contribution of this face to correction
                h    = face_center{i}(:,3);
                norm = face_center{i}(:,4);
                dF_int(j) = dF_int(j) + sum(pFace .* (dx*h) .* norm);
            end
        end
        
        % --------------------------------------------------------------------
        % Plot gap pressures
        % --------------------------------------------------------------------
        figure; 
        levels = -100:10:100;

        cmap_R1 = linspace(0.2, 0.7, 128);
        cmap_R2 = linspace(0.7, 0.6, 128);
        cmap_G1 = linspace(0.2, 0.7, 128);
        cmap_G2 = linspace(0.7, 0.0, 128);
        cmap_B1 = linspace(0.6, 0.7, 128);
        cmap_B2 = linspace(0.7, 0.2, 128);
        cmap = [cmap_R1', cmap_G1', cmap_B1'; cmap_R2', cmap_G2', cmap_B2'];

        subplot(3,1,1); box on; ylabel('y [m]'); hold on;
        [C1, h1] = contourf(zq, yq', p_int_LE{end-1}, levels);
        scatter(zLE_plot, yLE_plot, 10, 'w');
        title('x = 0 mm'); 
        caxis([levels(1) levels(end)]);
        colormap(cmap);
        set(gca,'xtick',[])

        subplot(3,1,2); box on; ylabel('y [m]'); hold on;
        [C2, h2] = contourf(zq, yq, p_int_TE{end-1}, levels);
        scatter(zTE_plot, yTE_plot, 10, 'w');
        title('x = 890 mm'); 
        caxis([levels(1) levels(end)]);
        colormap(cmap);
        set(gca,'xtick',[])

        subplot(3,1,3); box on; ylabel('y [m]'); xlabel('z [m]'); hold on;
        [C3, h3] = contourf(zq, yq, dp{end-1}, levels);
        title('\Delta p'); 
        caxis([levels(1) levels(end)]);
        colormap(cmap);

        clabel(C1,h1,levels);
        clabel(C2,h2,levels);
        clabel(C3,h3,levels);

    end

    % --------------------------------------------------------------------
    % Substitute interpolated data in output -----------------------------
    % --------------------------------------------------------------------

    P(:,[LE; TE]) = [];
    X([LE; TE],:) = [];

    P2 = zeros(size(P,1), length(zq)*length(yq));
    P3 = P2;
    X2 = zeros(length(zq)*length(yq),3);
    X3 = X2;

    for j=1:n1
        for k=1:n2

            for i=1:length(p_int_TE)
                P2(i,(j-1)*n2+k) = p_int_TE{i}(j,k);
                P3(i,(j-1)*n2+k) = p_int_LE{i}(j,k);
            end
            
            X2((j-1)*n2+k,:) = [ 445, yq(j), zq(k)];
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