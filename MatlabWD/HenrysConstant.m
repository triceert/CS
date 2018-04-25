function H = HenrysConstant(T,cmp)
% - Calculates Henry's constant for a given compound at a given temperature
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%        unt = unit struct
% OUTPUT: H: Henry's constant for the component [Pa]


%Assign Vectors from Compound struct
%   From Nitrogen to Hydrogen Cyanide (identifier 2-6)

kHstd = cmp(2:6).kHstd; %[M/atm]
deltaHsolR = cmp(2:6).deltaHsolR; %[K]
Tstd = 298.15; %[K]

kH = kHstd.*exp(-deltaHsolR.*(1/T-1/Tstd)); %[M/atm]
H = 55.3*101325/kH; %[Pa]
end


