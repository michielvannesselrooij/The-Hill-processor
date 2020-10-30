function wallCorr = hotwireWallCalibration(u, y)
%% Preparation

% Reshape data
x  = y;
y  = u;

% End values (x=0 and y=0)
y0 = y(end);
x1 = max(x);

% Initial conditions
a  = 0.0001;
b  = min(y)*1.5;
c  = 1.1;
d  = 0.96;

C0 = [a,    b,  c,   d];
lb = [1e-6, -1, 0.5, 0.5];
ub = [1,    0,  1.5, 1.5];

%% Optimization

% Function to be optimized (the near-wall weighted residual)
res = @(C) mean(abs( (C(1)*(y0 - C(2))./(x.^C(3) + C(1)) + C(2) - ...
    ( C(2)+C(1)*(y0-C(2))/(x1^C(3)+C(1)) )/( exp(x1^C(4))-1 ) * ...
    (exp(x.^C(4))-1) - y) .* 1./(x+1e-3) ));

% Optimization
options = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
C       = fmincon(res, C0, [], [], [], [], lb, ub, [], options);

%% Create output

a = C(1);
b = C(2);
c = C(3);
d = C(4);

A1 = a*(y0 - b);
A2 = -(b + A1/(x1^c+a))/(exp(x1^d)-1);

y1       = A1./(x.^c + a) + b;
y2       = A2*(exp(x.^d)-1);

wallCorr = y1+y2;

%% Plot results

figure; hold on; box on;
xlabel('y'); ylabel('wall correction [V]');
plot(x,y,'k:');
plot(x,wallCorr);
legend('data','fitted');

%%