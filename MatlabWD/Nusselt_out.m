function [Nu] = Nusselt_out(Re,Pr)
% INPUT: Pr Prandtl number
%        Re Reynolds number
% OUTPUT: Nu Nusselt number

Nu = 0.0288*Re^0.8*Pr^(1/3);
end

