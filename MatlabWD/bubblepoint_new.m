function [bubbleT] = bubblepoint_new(x, P, cmp, unt, thermo_model)
% calculates the bubble point for a binary mixture of H2O and HCN
% INPUT: x = vector of mole fractions; x(1) = H2O, x(2) = HCN; 
% P: pressure [Pa]
% cmp & unt: structs from Excel-file 


%x1 = x(1); 
options = optimset('Display', 'off'); 
bubbleT_0 = 50+273.15; 
delta_g12 = 500.9610; % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1)
delta_g21 = 539.9577; 
alpha12 = 0.3836; % this has nothing to do with the relative volatility alpha 
if strcmp(thermo_model, 'nrtl')
    gamma = @(bubbleT) nrtl(x(1), bubbleT, delta_g12, delta_g21, alpha12); 
elseif strcmp(thermo_model, 'vanlaar')
    gamma = @(bubbleT) vanlaar(x(1), bubbleT, delta_g12, delta_g21, alpha12); 
elseif strcmp(thermo_model, 'ideal')
    gamma = @(bubbleT) [1 1];   % ideal model, need to pass bubbleT to have same syntax as in other cases     
end
bubbleT = fsolve(@(bubbleT) bubbleT_solver(bubbleT, gamma,x,cmp), bubbleT_0, options); 





function [bT] = bubbleT_solver(bubbleT, gamma,x,cmp)
    gamma_vector = gamma(bubbleT);  
    gamma_HCN = gamma_vector(1); 
    gamma_H2O = gamma_vector(2); 
    P_sat = @(bubbleT) antoine_equation_new(cmp,bubbleT,n);
    bT = x(1).*P_sat(cmp,bubbleT,1).*gamma_H2O+x(2).*P_sat(cmp,bubbleT,6).*gamma_HCN-P; % species1 = HCN, species2 = H2O when it comes to gamma; 
end

end