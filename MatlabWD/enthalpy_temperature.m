function [HT] = enthalpy_temperature(T,cmp,unt)
% - Calculates the enthalpy of a GASEOUS compound with given parameters as a function of the temperature
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%        unt = unit struct
% OUTPUT: enthalpy [J.mol-1]

R=unt(5).idgc; % [J.mol-1.K-1]

boo = (T > 1000);

switch boo
    case 0 % T <= 1000K
        a = [cmp(1:6).a0_hcparlow;cmp(1:6).a1_hcparlow;cmp(1:6).a2_hcparlow;cmp(1:6).a3_hcparlow;cmp(1:6).a4_hcparlow];
        b1 = cmp(1:6).b1_hcparlow;
    case 1 % T > 1000K
        a = [cmp(1:6).a0_hcparhigh;cmp(1:6).a1_hcparhigh;cmp(1:6).a2_hcparhigh;cmp(1:6).a3_hcparhigh;cmp(1:6).a4_hcparhigh];
        b1 = cmp(1:6).b1_hcparhigh;
    otherwise
end

HT = R.*T.*(a(1,:) + a(2,:).*T/2 + a(3,:).*T^2/3 + a(4,:).*T^3/4 + a(5,:).*T^4/5)+R.*b1;
end

