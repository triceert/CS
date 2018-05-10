function [H] = HenrysConstant(T,cmp,n)
% - Calculates Henry's constant for a given compound at a given temperature
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%        n = compound number
% OUTPUT: H: Henry's constant for the component [Pa]


kHstd = cmp(n).kHstd; %[M/atm]
deltaHsolR = cmp(n).deltaHsolR; %[K]
Tstd = 298.15; %[K]

kH = kHstd*exp(-deltaHsolR*(1/T-1/Tstd)); %[M/atm]
H = 55.3*101325/kH; %[Pa]
end


