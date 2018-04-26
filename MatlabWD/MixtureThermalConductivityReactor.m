function [lambda_mix] = MixtureThermalConductivityReactor(F,T,plambda)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1:Nitrogen, 2:Methane, 3:Ammonia, 4: Hydrogen, 5:Hydrogen cyanide
%        T = temperature
%        plambda = 9x2 array with polynom of the compound on each line
%
% OUTPUT: lambda_mix [W.m-1.K-1]

lambdaT = plambda(2:6,1)*T + plambda(2:6,2);
F_tot = sum(F);
lambda_mix = sum(F.*lambdaT)/F_tot;
end