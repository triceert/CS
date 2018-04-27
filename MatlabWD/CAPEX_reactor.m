function [priceUSdollars] = CAPEX_reactor(number_of_tubes,unt)
% INPUT: number of tubes in the reactor    
%        unt = unit struct
% OUTPUT: price of the reactor in US$

price_tube = unt(1).price_tube;
priceUSdollars = price_tube*number_of_tubes;
end

