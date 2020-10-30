function E = muskerFit(s)
% CALCULATES RESIDUAL BETWEEN MUSKER PROFILE AND DESIRED B CONSTANT

% Load data
load('Musker_temp.mat', 'B', 'k', 'ypmax');

% Create Musker profile
yp_musker = 0:0.001:ypmax;
up_musker = zeros(size(yp_musker));
dy_musker = yp_musker(2)-yp_musker(1);
dudy      = ((yp_musker.^2)/k + 1/s) ./...
                                (yp_musker.^3 + (yp_musker.^2)/k + 1/s);
for i=2:length(yp_musker)
    up_musker(i) = up_musker(i-1) + dudy(i)*dy_musker;
end

% Determine constant B from Musker profile
B_musker = up_musker(end) - log(yp_musker(end))/k;

% Calculate error with desired value of B
E = abs(B-B_musker);