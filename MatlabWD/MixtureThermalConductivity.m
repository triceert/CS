function [lambda_mix] = MixtureThermalConductivity(F,T,plambda)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1: Water, 2:Nitrogen,3:Methane, 4:Ammonia,5: Hydrogen, 6:Hydrogen cyanide,
% 7:Sulfuric acid, 8:Ammonium sulfate, 9:Carbon dioxide
%        T = temperature
%        plambda = 9x2 array with polynom of the compound on each line
%
% OUTPUT: lambda_mix [W.m-1.K-1]

lambdaT = plambda(:,1)*T + plambda(:,2);
F_tot = sum(F);
lambda_mix = sum(F.*lambdaT)/F_tot;
end