function DR = logShiftToDrag(dB, k, rho, V, u_tau)
% ESTIMATE THE DRAG EFFECT OF A SHIFT IN THE LOG-LAYER CONSTANT B
% Based on Tani 1988, specifically equation 4 and accompanying text
% ---------------------------------------
%
% INPUTS
% ---
% dB        : shift in log-layer constant B     Matrix of size [M,N]
% k         : von Karman constant               Scalar or matrix of size [M,N]
% rho       : flow density [kg/m3]              Scalar or matrix of size [M,N]
% V         : freestream velocity [m/s]         Scalar or matrix of size [M,N]
% u_tau     : friction velocity [m/s]           Scalar or matrix of size [M,N]
%
% OUTPUTS
% ---
% DR        : drag reduction in percents        Matrix of size [M,N]
%
% ---------------------------------------

%% Check inputs

if size(k) ~= size(dB)
    if length(k) ~= 1
        error('k should be a scalar or equal to the size of dB');
    end
end

if size(rho) ~= size(dB)
    if length(rho) ~= 1
        error('rho should be a scalar or equal to the size of dB');
    end
end

if size(V) ~= size(dB)
    if length(V) ~= 1
        error('V should be a scalar or equal to the size of dB');
    end
end

if size(u_tau) ~= size(dB)
    if length(u_tau) ~= 1
        error('u_tau should be a scalar or equal to the size of dB');
    end
end

%% Calculation

tau_w  = u_tau.^2 * rho;
cf     = tau_w ./ (1/2*rho.*V.^2);

DR     = 2*k*dB./(1+k*sqrt(2./cf)) * 100;