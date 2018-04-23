function [P_sat] = antoine_equation(A, B, C, T)
% - Calculates saturation pressure of a compound with given Antoine
%   parameters for a given temperature
% - Can also submit vectors of A, B, C (T can be scalar --> all components
%   at same temperature) to find the saturation pressures for the individual
%   components (not yet calculated for mixture)
% INPUT: A, B, C = Antoine parameters (mmHg units)
%        T = Temperature [K]
P_sat_unconverted = exp(A-B./(T+C));    % saturation pressure [mmHg]
P_sat = P_sat_unconverted.*133.3;       % 1 mmHg = 133.3 Pa
end

