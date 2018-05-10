function [unt] = TOTEX_reactor(unt)
% INPUT: unt = unit struct
% OUTPUT: price of the reactor in US$

CAPEX = unt(1).capex; %[$]
OPEX = unt(1).opex; %[$.y-1]
TOTEX = CAPEX + 10*OPEX;  %for 10 years

%CAPEX in unt
unt(1).totex = TOTEX;
end

