function E = canonicalFit(x)
% CALCULATES THE RESIDUAL OF THE FIT OF MEAN BOUNDARY LAYER DATA WITH A
% CANONICAL MODEL BASED ON THE APPROACH BY LOPEZ, 2015.

if length(x) ~= 5
    error(['canonicalFit expects a vector of length 5: u_tau, dy, k,'...
        'PI, and d. Check your input.']);
end

% ---------------------
% Unpack inputs
% ---------------------
ut = x(1);
dy = x(2);
k  = x(3);
PI = x(4);
d  = x(5);

% ---------------------
% Load experiment data
% ---------------------
load('BL_temp.mat', 'y', 'u', 'nu');

% ---------------------
% Create model
% ---------------------
[yp_model, up_model] = canonicalBL(ut, k, PI, d, nu);

% -------------------
% Determine error
% -------------------
up = u/ut;
yp = (y+dy)*ut/nu;
up_model_int = interp1(yp_model, up_model, yp, 'linear', 'extrap');

E = mean(abs(up_model_int - up)./up_model_int);
% E = mean(sqrt((up_model_int - up).^2));