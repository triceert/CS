function [lambda_mix] = MixtureThermalConductivityNaturalGas(x,T,plambda)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: x = molar fraction of gas : 1:Water, 2:CO2
%        T = temperature
%        plambda = 9x2 array with polynom of the compound on each line
%
% OUTPUT: lambda_mix [W.m-1.K-1]

lambdaT = plambda([1 9],1)*T + plambda([1 9],2);
lambda_mix = sum(x.*lambdaT);
end