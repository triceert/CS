function [mu_mix] = MixtureDynamicViscosityNaturalGas(x,T,pmu)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: x = molar fraction of gas : 1:Water, 2:CO2
%        T = temperature
%        pmu = 9x3 array with polynom of the compound on each line [Pa.s]
%
% OUTPUT: mu_mix

muT = pmu([1 9],1)*T^2 + pmu([1 9],2)*T + pmu([1 9],3);
mu_mix = sum(x.*muT);
end