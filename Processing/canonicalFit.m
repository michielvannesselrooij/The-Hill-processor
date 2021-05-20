function E = canonicalFit(x)
% CALCULATES THE RESIDUAL OF THE FIT OF MEAN BOUNDARY LAYER DATA WITH A
% CANONICAL MODEL BASED ON THE APPROACH BY LOPEZ, 2015.

% -----------------------
% Load normalization data
% -----------------------
if exist('norm.mat','file') == 2
    
    load('norm.mat', 'x0', 'E0');
    
else % If does not exist, do not normalize (scale by 1)
    
    x0 = ones(5,1);
    E0 = 1;
    
end

% -----------------------
% Unpack inputs
% -----------------------
if length(x) == 5
    
    ut = x(1) * x0(1);
    dy = x(2) * x0(2);
    k  = x(3) * x0(3);
    PI = x(4) * x0(4);
    d  = x(5) * x0(5);
        
elseif isstruct(x) 
    
    ut = x.ut * x0(1);
    dy = x.dy * x0(2);
    k  = x.k  * x0(3);
    PI = x.PI * x0(4);
    d  = x.d  * x0(5);
        
else
    
    error(['canonicalFit expects a vector of length 5: u_tau, dy, k,'...
        'PI, and d. Check your input.']);
    
end

% -----------------------
% Load experiment data
% -----------------------
load('BL_temp.mat', 'y', 'u', 'nu', 'muskerOpt');

% -----------------------
% Create model
% -----------------------
[yp_model, up_model] = canonicalBL(ut, k, PI, d, nu, muskerOpt);

% ---------------------
% Determine error
% ---------------------
up           = u/ut;
yp           = (y+dy)*ut/nu;
up_model_int = interp1(yp_model, up_model, yp, 'linear', 'extrap');

% Check and penalize potential negative values
up_model_int(up_model_int <= 0) = 1e-3;

E = mean(abs(up_model_int - up)./up_model_int);
% E = mean(sqrt((up_model_int - up).^2));

% Normalize error
E = E/E0;