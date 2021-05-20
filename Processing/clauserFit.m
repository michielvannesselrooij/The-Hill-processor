function E = clauserFit(xn)
% Optimization function for Clauser fit of log-layer
% -----------------------------------------------------

% Load temporary data
load('BL_clauser_temp.mat', 'y', 'u', 'nu');

if exist('BL_clauser_temp_norm.mat','file') == 2
    
    load('BL_clauser_temp_norm.mat', 'x0', 'E0');
    
else % If does not exist, do not normalize (scale by 1)
    
    x0 = ones(3,1);
    E0 = 1;
    
end

% Unpack inputs
x  = xn .* x0;
k  = x(1);
B  = x(2);
ut = x(3);

% Non-dimensionalize profile
yp  = y*ut/nu;
up  = u/ut;

% Select y range
yp_min = 100;
yp_max = 500;

[idx, ~] = find(yp>yp_min & yp<yp_max);

yp  = yp(idx);
up  = up(idx);

% Stop if too few data points
if length(idx) < 5
    E = 2*E0;
    return;
end

% Calculate Clauser profile
up_Clauser = 1/k * log(yp) + B;

% Determine normalized error
E = mean(abs(up_Clauser - up));
E = E/E0;