function [ut, dy, k, B, PI, d, yp_model_q, up_model_q] = ...
                                        canonicalOpt(y, u, nu, showPlot)
% CALCULATES THE PROPERTIES OF MEAN BOUNDARY LAYER DATA WITH A
% CANONICAL MODEL BASED ON THE APPROACH BY LOPEZ, 2015.

% ---------------------------------------------------------------------
% Initial conditions
% ---------------------------------------------------------------------

% Ut
k   = 0.41;
B   = 5.0;
U   = u(end);
y(y==0) = 1e-9; % Prevent ln(0)=-Inf
res = @(Cf) mean(abs( ((1/k*sqrt(Cf/2))*log(y*U/nu) + ...
                1/k*sqrt(Cf/2)*log(sqrt(Cf/2)) + B*sqrt(Cf/2)) - u/U ));

options    = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
Cf         = fmincon(res, 1, [], [], [], [], 1e-9, 10, [], options);
ut_clauser = sqrt(Cf/2)*U;

ut_min = 0.01;
ut_max = 2*ut_clauser;
ut_0   = ut_clauser;

% dy
dy_min = -1e-3;
dy_max = 5e-3;
dy_0   = 0;

%k
k_min  = 0.1;
k_max  = 0.9;
k_0    = 0.4;

% PI
PI_min = 0;
PI_max = 5;
PI_0   = 0.5;

% d
d99   = interp1(u/u(end), y, .99);
d_min = 0.1*d99;
d_max = 2.0*d99;
d_0   = 1.3*d99;

% ---------------------------------------------------------------------
% Optimization
% ---------------------------------------------------------------------

% Temporarily store boundary layer data for optimization function to use
save('BL_temp.mat', 'y', 'u', 'nu');

% Fit optimization parameters
x0 = [ut_0,   dy_0,   k_0,   PI_0,   d_0];
lb = [ut_min, dy_min, k_min, PI_min, d_min];
ub = [ut_max, dy_max, k_max, PI_max, d_max];

% Estimated value of s parameter
options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp',...
            'MaxFunctionEvaluations',1000, 'FunctionTolerance', 1e-10);
x = fmincon(@canonicalFit, x0, [], [], [], [], lb, ub, [], options);

% Optimized value of s parameter
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp',...
%             'MaxFunctionEvaluations',1000);
% lb = x*0.9; ub = x*1.1;
% x = fmincon(@canonicalFit2, x, [], [], [], [], lb, ub, [], options);

ut = x(1);
dy = x(2);
k  = x(3);
PI = x(4);
d  = x(5);

% Remove temporary data file
delete('BL_temp.mat');

% Resample model BL at reasonable size
[yp_model, up_model, B] = canonicalBL(ut, k, PI, d, nu);
yp_model_q = logspace(0,4,1000);
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