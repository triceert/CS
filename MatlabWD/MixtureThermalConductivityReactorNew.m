function [lambda_mix] = MixtureThermalConductivityReactorNew(cmp,F,T)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1:Nitrogen, 2:Methane, 3:Ammonia, 4: Hydrogen, 5:Hydrogen cyanide
%        T = temperature
%
% OUTPUT: lambda_mix [W.m-1.K-1]

A = extractfield(cmp(2:6),'condA')';
B = extractfield(cmp(2:6),'condB')';
C = extractfield(cmp(2:6),'condC')';

lambdaT = (C*T^2 + B*T + A);
F_tot = sum(F);
lambda_mix = sum(F.*lambdaT)/F_tot;
end