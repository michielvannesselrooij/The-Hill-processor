function [Re, V, rho, nu] = velocityCalculator(q, pa, T, hum)
% ------------------------------------------------------------------------
% [For velocity sweep] use pitot pressure, absolute pressure, temperature
% and humidity to calculate density, velocity and Reynolds number 
% 
% MvN 2019 - Dimple Aerospace BV
% ------------------------------------------------------------------------
    
    % Prepare
    q      = q-q(1);
    q(q<0) = NaN;

    % Calculate velocities
    [Re, V, rho, nu] = calcV(q, pa, T, hum);

end