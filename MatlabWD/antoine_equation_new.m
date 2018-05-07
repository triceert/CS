function [P_sat] = antoine_equation_new(cmp,T,n)
% Calculates saturation pressure of a compound with given Antoine
%   parameters for a given temperature
% INPUT: cmp = compound struct
%        T = temperature [K]
%        n = compound index
% OUTPUT: saturation pressure [Pa]

boo1 = (T > 373);

switch boo1
    case 0 %(T <= 100Cels)
        A = extractfield(cmp(1:7),'antaRT')';
        B = extractfield(cmp(1:7),'antbRT')';
        C = extractfield(cmp(1:7),'antcRT')';
    case 1 %(T > 100Cels)
        A = extractfield(cmp(1:7),'anta100')';
        B = extractfield(cmp(1:7),'antb100')';
        C = extractfield(cmp(1:7),'antc100')';
end

boo2 = (n == 7); %Different source for sulfuric acid

switch boo2
    case 0 %(n <= 6)
        T = T -273.15; % Formula in �C, but all T always handled in K, so need to convert
        P_sat_unconverted = 10.^(A-B./(T+C)); %in psi (pounds per square inch)
        P_sat = P_sat_unconverted * 6894.76; 
    case 1 %(n = 7)
       P_sat = exp(A-B./(T+C)); 
end

end

