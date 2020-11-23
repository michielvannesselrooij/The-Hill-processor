function [yp_model, up_model, B] = canonicalBL2(ut, k, PI, d, nu)
%CREATES CANONICAL MODEL BASED ON THE APPROACH BY LOPEZ, 2015.

% Calculate integration constant s
B  = 5.0;
lb = -10;
ub = 10;
r  = @(B) 1.6*(exp(0.1663*B)-1) - k*B;
options = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');        
B = fmincon(r, B, [], [], [], [], lb, ub, [], options);

% Prepare variables
ypmax = d*ut/nu;           % y_plus @ boundary layer edge
yp1   = (0:0.001:ypmax)';  % y_plus in boundary layer

% Part 1: Musker profile

    % Option 1: fix for k=0.41, B=5.0
    % s = 0.001093;

    % First estimate s (MvN method)
    C_a = [-0.0046, -0.4951, 1.0933, -0.0509];
    C_b = [0.1343,  -0.6548, 1.0593, -0.7775, 0.6402];
    C_c = [-0.7951, 1.4998,  0.8492];
    a = polyval(C_a,k);
    b = polyval(C_b,k);
    c = C_c(1)./k.^C_c(2) + C_c(3);

    s = (a/(B-c))^(1/b);

    % Optimize fit with B
    save('Musker_temp.mat', 'B', 'k', 'ypmax');
    options = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
    s = fmincon(@muskerFit, s, [], [], [], [], [], [], [], options);
    delete('Musker_temp.mat');

    yp_musker = yp1;
    up_musker = zeros(size(yp_musker));
    dy_musker = yp_musker(2)-yp_musker(1);
    dudy      = ((yp_musker.^2)/k + 1/s) ./...
                                    (yp_musker.^3 + (yp_musker.^2)/k + 1/s);
    for i=2:length(yp_musker)
        up_musker(i) = up_musker(i-1) + dudy(i)*dy_musker;
    end
    B = up_musker(end) - log(yp_musker(end))/k;

% Part 2: Bump profile
yp_bump = yp1;
M1      = 30;
M2      = 2.85;
up_bump = exp(-log(yp_bump/M1).^2)/M2;

% Part 3: Outer profile
yp_outer = yp1;
eta      = 0.001:0.001:1;
a2       = 132.8410;
a3       = -166.2041;
a4       = 71.9114;
W_outer  = ((1-exp(-0.25*(5*a2+6*a3+7*a4)*eta.^4 + a2*eta.^5 + ...
            a3*eta.^6 + a4*eta.^7))./ (1-exp(-0.25*(a2 + 2*a3 + 3*a4)))) ...
            .* (2*PI-log(eta))/k;
W_outer  = interp1(eta*d*ut/nu, W_outer, yp_outer, 'linear', 'extrap'); 

% Combine model elements
up1      = up_musker + up_bump + W_outer;
up_model = [up1; up1(end)];  % Extend beyond boundary layer
yp_model = [yp1; 1e6*ypmax]; % Extend beyond boundary layer
