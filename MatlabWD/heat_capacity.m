function [cpT] = heat_capacity(T,cmp,unt)
% - Calculates the heat capacity of a GASEOUS compound with given
% parameters as a function of the temperature
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%        unt = unit struct
% OUTPUT: heat capacity [J.mol-1.K-1]

R=unt(5).idgc; % [J.mol-1.K-1]

boo = (T > 1000);

switch boo
    case 0 % T <= 1000K
        a = [cmp(1:6).a0_hcparlow;cmp(1:6).a1_hcparlow;cmp(1:6).a2_hcparlow;cmp(1:6).a3_hcparlow;cmp(1:6).a4_hcparlow];
    case 1 % T > 1000K
        a = [cmp(1:6).a0_hcparhigh;cmp(1:6).a1_hcparhigh;cmp(1:6).a2_hcparhigh;cmp(1:6).a3_hcparhigh;cmp(1:6).a4_hcparhigh];
    otherwise
end


cpT = R*(a(1,:) + a(2,:).*T + a(3,:).*T.^2 + a(4,:).*T.^3 + a(5,:).*T.^4);

end