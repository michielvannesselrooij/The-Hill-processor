function [Re, V, rho, nu] = velocityCalculator(q, pa, T, hum)
% ------------------------------------------------------------------------
% Uses pitot pressure, absolute pressure, temperature and humidity to
% calculate density, velocity and Reynolds number
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------
    
    % Prepare
    T0     = 273.15;
    q      = q-q(1);
    q(q<0) = NaN;

    % Vapor pressure following Arden Buck 1996
    p_vap = hum .*0.61121 .*exp( (18.678-T/234.84).*(T./(257.14+T)) ) /1000;
    
    % Calculate properties of interest
    rho = ( (pa-p_vap)*0.028964 + p_vap*0.018016 )./(8.314.*(T+T0));
    nu  = 1./pa *4.18528*1e-4 .*(T+T0).^2.5 ./(110.4+T+T0);
    V   = sqrt(2*q./rho);
    Re  = 1*V./nu;

end