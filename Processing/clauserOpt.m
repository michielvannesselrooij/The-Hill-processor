function [k, B, ut] = clauserOpt(y, u, nu, k0, B0, ut0, showPlot)
% Fit boundary layer data to log layer region
% -----------------------------------------------------

% Set and scale initial conditions and bounds
x0 = [k0,  B0, ut0];
lb = [0.1, -5, 0.01];
ub = [1.0, 15, 15];

x0n  = x0./x0;
lbn  = lb./x0;
ubn  = ub./x0;

% Temporarily store data for optimizer
save('BL_clauser_temp.mat', 'y', 'u', 'nu');

% Initial error
E0 = clauserFit(x0);

% Temporarily store normalization data
save('BL_clauser_temp_norm.mat', 'x0', 'E0');

% Run optimizer
disp('Running Clauser profile fitter...');
options = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
xn = fmincon(@clauserFit, x0n, [], [], [], [], lbn, ubn, [], options);
disp(' ');

% Unpack results
x  = xn .* x0;
k  = x(1);
B  = x(2);
ut = x(3);

% Clean up
delete('BL_clauser_temp.mat');
delete('BL_clauser_temp_norm.mat');

% Show plot (optional)
if showPlot
    
    yp         = y*ut/nu;
    up         = u/ut;

    yp_min     = 100;
    yp_max     = 500;

    [idx, ~]   = find(yp>yp_min & yp<yp_max);

    yp_Clauser = yp(idx);
    up_Clauser = 1/k * log(yp_Clauser) + B;
    
    figure; hold on;
    plot(yp, up, 'ks-');
    plot(yp_Clauser, up_Clauser, 'r-','LineWidth',1.5)
    legend('Data', 'Clauser fit','Location','NorthWest');
    xlabel('y^+');
    ylabel('u^+');
    set(gca, 'XScale', 'log')
    
end