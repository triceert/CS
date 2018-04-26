function [mu_mix] = MixtureDynamicViscosityReactor(F,T,pmu)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1:Nitrogen, 2:Methane, 3:Ammonia, 4: Hydrogen, 5:Hydrogen cyanide
%        T = temperature
%        pmu = 9x3 array with polynom of the compound on each line [Pa.s]
%
% OUTPUT: mu_mix

muT = pmu(2:6,1)*T^2 + pmu(2:6,2)*T + pmu(2:6,3);
F_tot = sum(F);
mu_mix = sum(F.*muT)/F_tot;
end