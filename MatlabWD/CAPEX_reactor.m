function [CAPEX_reactor_priceUSdollars] = CAPEX_reactor(unt)
% INPUT: unt = unit struct
% OUTPUT: price of the reactor in US$

price_tube = unt(1).price_tube;
number_of_tubes = unt(1).N_tubes;
CAPEX_reactor_priceUSdollars = price_tube*number_of_tubes;
end

