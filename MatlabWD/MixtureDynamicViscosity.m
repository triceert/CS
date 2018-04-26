function [mu_mix] = MixtureDynamicViscosity(F,T,pmu)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1: Water, 2:Nitrogen,3:Methane, 4:Ammonia,5: Hydrogen, 6:Hydrogen cyanide,
% 7:Sulfuric acid, 8:Ammonium sulfate, 9:Carbon dioxide
%        T = temperature
%        pmu = 9x3 array with polynom of the compound on each line [Pa.s]
%
% OUTPUT: mu_mix

muT = pmu(:,1)*T^2 + pmu(:,2)*T + pmu(:,3);
F_tot = sum(F);
mu_mix = sum(F.*muT)/F_tot;
end