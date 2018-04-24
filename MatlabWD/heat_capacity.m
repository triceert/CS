function [cpT] = heat_capacity(T,a)
% - Calculates the heat capacity of a GASEOUS compound with given 
%
% INPUT: a = coefficients from the NASA as 1x5 or 5x1 vector
%        T = Temperature [K]

R = 8.314; % [J.mol-1.K-1]
cpT = R*(a(1) + a(2)*T + a(3)*T^2 + a(4)*T^3 + a(5)*T^4); % [J.mol-1.K-1]

end