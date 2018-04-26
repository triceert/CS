function [cpT] = heat_capacity(T,cmp,unt,n)
% - Calculates the heat capacity of a GASEOUS compound with given
% parameters as a function of the temperature
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%        unt = unit struct
%        n idientifier of compound in excel table
%        1   {'H2O';'N2';'CH4';'NH3';'H2';'HCN'} 6
% OUTPUT: heat capacity [J.mol-1.K-1]

R=unt(5).idgc; % [J.mol-1.K-1]

boo = (T > 1000);

switch boo
    case 0 % T <= 1000K
        a = [cmp(n).a0_hcparlow;cmp(n).a1_hcparlow;cmp(n).a2_hcparlow;cmp(n).a3_hcparlow;cmp(n).a4_hcparlow];
    case 1 % T > 1000K
        a = [cmp(n).a0_hcparhigh;cmp(n).a1_hcparhigh;cmp(n).a2_hcparhigh;cmp(n).a3_hcparhigh;cmp(n).a4_hcparhigh];
    otherwise
end


cpT = R*(a(1,:) + a(2,:).*T + a(3,:).*T.^2 + a(4,:).*T.^3 + a(5,:).*T.^4);



end