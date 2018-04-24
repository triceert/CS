function [HT] = enthalpy_temperature(T,a,b1)
% - Calculates the enthalpy of a GASEOUS compound with given parameters as a function of the temperature
%
% INPUT: a = coefficients a0,a1,a2,a3,a4 from the NASA as 1x5 or 5x1 vector
%        b1 = coefficient b1 from the NASA
%        T = Temperature [K]

R = 8.314; % [J.mol-1.K-1]
HT = RT*(a(1) + a(2)*T/2 + a(3)*T^2/3 + a(4)*T^3/4 + a(5)*T^4/5)+R*b1; % [J.mol-1]
end

