function [Nu] = Nusselt_in(Re,Pr)
% INPUT: Pr Prandtl number
%        Re Reynolds number
% OUTPUT: Nu Nusselt number

Nu = 0.023*Re^0.8*Pr^0.4;
end

