function [ut, dy, k, B, PI, d, yp_model_q, up_model_q] = ...
                                        canonicalOpt(y, u, nu, showPlot)
% CALCULATES THE PROPERTIES OF MEAN BOUNDARY LAYER DATA WITH A
% CANONICAL MODEL BASED ON THE APPROACH BY LOPEZ, 2015.

% ---------------------------------------------------------------------
% Initial conditions
% ---------------------------------------------------------------------

% Ut
k          = 0.41;
B          = 5.0;
U          = u(end);
y(y==0)    = 1e-9; % Prevent ln(0)=-Inf

res        = @(Cf) mean(abs( ( (1/k*sqrt(Cf/2))*log(y*U/nu) + ...
                               (1/k*sqrt(Cf/2))*log(sqrt(Cf/2)) + B*sqrt(Cf/2) ) - u/U ));

options    = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
Cf         = fmincon(res, 1, [], [], [], [], 1e-8, 10, [], options);
ut_Clauser = sqrt(Cf/2)*U;

% Estimate ut with Clauser method
% k0         = 0.41;
% B0         = 5.0;
% ut0        = 1;
% showPlot   = 0;
% [~, ~, ut_Clauser] = clauserOpt(y, u, nu, k0, B0, ut0, showPlot);

ut_min     = 0.01;
ut_max     = 1.5*ut_Clauser;
ut_0       = 1.3*ut_Clauser;

% dy
dy_min     = -min(y) + 1e-9;
dy_max     = 5e-3;
dy_0       = 1e-5; % prevent division by zero

% k
k_min      = 0.2;
k_max      = 0.8;
k_0        = 0.4;

% PI
PI_min     = 0;
PI_max     = 5;
PI_0       = 0.5;

% d
d99        = interp1(u/u(end), y, .99);
d_min      = 0.5*d99;
d_max      = 8.0*d99;
d_0        = 1.0*d99;

% ---------------------------------------------------------------------
% Optimization
% ---------------------------------------------------------------------

% Temporarily store boundary layer data for optimization function to use
muskerOpt = false;
save('BL_temp.mat', 'y', 'u', 'nu', 'muskerOpt');

x0 = [ut_0,   dy_0,   k_0,   PI_0,   d_0];
lb = [ut_min, dy_min, k_min, PI_min, d_min];
ub = [ut_max, dy_max, k_max, PI_max, d_max];

% Normalize inputs
x0n  = x0./x0;
lbn  = lb./x0;
ubn  = ub./x0;

% Find initial error
E0   = canonicalFit(x0);

% Temperarily store normalization data
save('norm.mat', 'x0', 'E0');

% Run optimizer
disp('Running Rodriguez-Lopez optimizer...');
options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp',...
            'MaxFunctionEvaluations', 1000, 'FunctionTolerance', 1e-9, ...
            'StepTolerance', 1e-9);
xn = fmincon(@canonicalFit, x0n, [], [], [], [], lbn, ubn, [], options);
disp(' ');

% Unpack results
x  = xn .* x0;
ut = x(1);
dy = x(2);
k  = x(3);
PI = x(4);
d  = x(5);

% Remove temporary data files
delete('BL_temp.mat');
delete('norm.mat');

% Resample model BL at reasonable size
muskerOpt = 1;
[yp_model, up_model, B] = canonicalBL(ut, k, PI, d, nu, muskerOpt);
yp_model_q = logspace(0, 4, 500);
up_model_q = interp1(yp_model, up_model, yp_model_q, 'linear', 'extrap');

% -------------------
% Plot result
% -------------------
if exist('showPlot','var')
    if showPlot==1
        yp = y*ut/nu;
        up = u/ut;

        figure; hold on; grid on; box on;
        h1 = plot(yp_model_q, up_model_q, 'k-');
        h2 = plot(yp, up, 'ko');
        
        set(gca,'xscale','log');
        xlim([1 2*max(yp)]);
        ylim([0 1.1*max(up)]);
        xlabel('y^+');
        ylabel('u^+');
        legend([h1, h2], 'model', 'experiment')
    end
end