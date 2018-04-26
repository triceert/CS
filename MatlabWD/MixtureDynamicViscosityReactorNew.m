function [mu_mix] = MixtureDynamicViscosityReactorNew(cmp,F,T)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: F = flux with indexes : 1:Nitrogen, 2:Methane, 3:Ammonia, 4: Hydrogen, 5:Hydrogen cyanide
%        T = temperature
%
% OUTPUT: mu_mix

A = extractfield(cmp(2:6),'viscA')';
B = extractfield(cmp(2:6),'viscB')';
C = extractfield(cmp(2:6),'viscC')';

muT = (C*T^2 + B*T + A)'*10^(-6);
F_tot = sum(F);
mu_mix = sum(F.*muT)/F_tot;
end