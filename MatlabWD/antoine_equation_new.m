function [P_sat] = antoine_equation_new(cmp,T,n)
% Calculates saturation pressure of a compound with given Antoine
%   parameters for a given temperature
% INPUT: cmp = compound struct
%        T = temperature [K]
%        n = compound index
% OUTPUT: P_sat = saturation pressure [Pa]

A = extractfield(cmp(n),'anta')';
B = extractfield(cmp(n),'antb')';
C = extractfield(cmp(n),'antc')';

t = T - 273.15; % Formula in Celsius, but all T always handled in K, so need to convert
P_sat_unconverted = 10.^(A-B./(t+C)); %in mmHg (pounds per square inch)
P_sat = P_sat_unconverted * 133.322; %%% convert to Pa

end

