function [d, d_star, theta, H] = calcBLproperties(u, y, U, d)
% ---------------------------------------------------------------------
% Calculates basic boundary layer properties
%
% INPUT:
% (required) u - velocity vector
% (required) y - position vector
% (optional) U - freestream velocity
% (optional) d - boundary layer thickness
% 
% MvN 2019 - Dimple Aerospace BV
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Error checking
% ---------------------------------------------------------------------
if iscell(u) || iscell(y)
    error('Input arguments u and y should not be cells');
end

if min(size(u)) ~= 1 || min(size(y)) ~= 1 || max(size(u)) ~= max(size(y))
    error(['Input arguments u and y should be 1-dimensional '...
        'matrices of equal length']);
end

% ---------------------------------------------------------------------
% Prepare data
% ---------------------------------------------------------------------
if size(y,1) > 1
    y = y';
end

if size(u,1) > 1
    u = u';
end

if y(1) > 0
    y = [0 y];
    u = [0 u];
end

if ~exist('U','var')
    U = max(u);
end

% ---------------------------------------------------------------------
% d:      Boundary layer thickness
% ---------------------------------------------------------------------
y_d    = .99;
if ~exist('d','var')
    d      = interp1(u/u(end), y, y_d);
end

% ---------------------------------------------------------------------
% Interpolation of data
% ---------------------------------------------------------------------
N      = 10000;
y_int  = d*(1/N : 1/N : y_d);
u_int  = interp1(y, u, y_int, 'linear', 'extrap');

% ---------------------------------------------------------------------
% d_star: Displacement thickness
% ---------------------------------------------------------------------
d_star = sum( (1 - u_int/U)*d/N );

% ---------------------------------------------------------------------
% theta:  Momentum thickness
% ---------------------------------------------------------------------
theta  = sum( (u_int/U .* (1 - u_int/U))*d/N );

% ---------------------------------------------------------------------
% H:      Shape factor
% ---------------------------------------------------------------------
H      = d_star/theta;
