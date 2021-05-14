%% USE AVAILABLE u AND y DATA

% Conditions
% ut = [1.15 1.1];
nu = 1.51e-5;

% Viscous sublayer model
visc_y = 1:0.1:10;
visc_u = visc_y;

% Log-layer model
B = 3;
k = 0.35;
log_y = 10:10:500;
log_u = 1/k * log(log_y) + B;

% U_tau with Clauser method
for i=1:length(u)
    for j=1:length(u{i})
        k               = 0.41;
        B               = 5.0;
        U               = u{i}{j}(end);
        y{i}{j}(y{i}{j}==0)   = 1e-9; % Prevent ln(0)=-Inf
        res             = @(Cf) mean(abs( ((1/k*sqrt(Cf/2))*log(y{i}{j}*U/nu) + ...
                          1/k*sqrt(Cf/2)*log(sqrt(Cf/2)) + B*sqrt(Cf/2)) - u{i}{j}/U ));

        options         = optimoptions(@fmincon, 'Display', 'none', 'Algorithm', 'sqp');
        Cf              = fmincon(res, 1, [], [], [], [], 1e-9, 10, [], options);
        ut(i,j)         = sqrt(Cf/2)*U
    end
end

for i=1:length(u)
    
    % Prepare plot
    figure; hold on; box on; grid on;
    set(gca,'xscale','log');
    c='kr';
    
    % Plot model parts
    plot(visc_y, visc_u, 'k-');
    plot(log_y, log_u, 'k-');
    
    % Plot profiles
    for j=1:length(u{i})
        up{i}{j} = u{i}{j}/ut(j);
        yp{i}{j} = y{i}{j}*ut(j)/nu;
        
        h(j)=plot(yp{i}{j},up{i}{j},'s-','Color',c(j));
    end
    
    legend(h,'Ref.','Target');
end


