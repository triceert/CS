function [P_sat] = antoine_equation_new(cmp,T,n)
% Calculates saturation pressure of a compound with given Antoine
%   parameters for a given temperature
% INPUT: cmp = compound struct
%        T = temperature [K]
%        n = compound index
% OUTPUT: P_sat = saturation pressure [Pa]

boo1 = (T > 373);

switch boo1
    case 0 %(T <= 100Cels)
        A = extractfield(cmp(n),'antaRT')';
        B = extractfield(cmp(n),'antbRT')';
        C = extractfield(cmp(n),'antcRT')';
    case 1 %(T > 100Cels)
        A = extractfield(cmp(n),'anta100')';
        B = extractfield(cmp(n),'antb100')';
        C = extractfield(cmp(n),'antc100')';
end

boo2 = (n == 7); %Different source for sulfuric acid

switch boo2
    case 0 %(n <= 6)
        T = T -273.15; % Formula in Celsius, but all T always handled in K, so need to convert
        P_sat_unconverted = 10.^(A-B./(T+C)); %in psi (pounds per square inch)
        P_sat = P_sat_unconverted * 6894.76; %%% double check this
    case 1 %(n = 7)
       P_sat = exp(A-B./(T+C)); 
end

end

