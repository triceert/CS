function [lambda_mix] = MixtureThermalConductivityNaturalGasNew(cmp,x,T)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: x = molar fraction of gas : 1:Water, 2:CO2
%        T = temperature
%       
% OUTPUT: lambda_mix [W.m-1.K-1]

A = [cmp(1).condA;cmp(9).condA];
B = [cmp(1).condB;cmp(9).condB];
C = [cmp(1).condC;cmp(9).condC];

lambdaT = (C*T^2 + B*T + A);
lambda_mix = sum(x.*lambdaT);
end