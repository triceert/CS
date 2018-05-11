function [mu_mix] = MixtureDynamicViscosityNaturalGasNew(cmp,x,T)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: x = molar fraction of gas : 1:Water, 2:CO2
%        T = temperature
%        
% OUTPUT: mu_mix [Pa.s]

A = [cmp(1).viscA;cmp(9).viscA];
B = [cmp(1).viscB;cmp(9).viscB];
C = [cmp(1).viscC;cmp(9).viscC];

muT = (C*T^2 + B*T + A)*10^(-7);
mu_mix = sum(x.*muT);
end