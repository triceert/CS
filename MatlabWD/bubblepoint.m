function [bubbleT] = bubblepoint(x, P, cmp, unt, thermo_model)
% calculates the bubble point for a binary mixture of H2O and HCN
% INPUT: x = vector of mole fractions; x(1) = H2O, x(2) = HCN; 
% P: pressure [Pa]
% cmp & unt: structs from Excel-file 


x1 = x(1); 
options = optimset('Display', 'off'); 
bubbleT_0 = 50+273.15; 
delta_g12 = 500.9610; % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1)
delta_g21 = 539.9577; 
alpha12 = 0.3836; % this has nothing to do with the relative volatility alpha 
if strcmp(thermo_model, 'nrtl')
    gamma = @(bubbleT) nrtl(x1, bubbleT, delta_g12, delta_g21, alpha12); 
elseif strcmp(thermo_model, 'vanlaar')
    gamma = @(bubbleT) vanlaar(x1, bubbleT, delta_g12, delta_g21, alpha12); 
end
bubbleT = fsolve(@(bubbleT) bubbleT_solver(bubbleT, gamma), bubbleT_0, options); 






function [bT] = bubbleT_solver(bubbleT, gamma)
    gamma_vector = gamma(bubbleT);  
    gamma_HCN = gamma_vector(1); 
    gamma_H2O = gamma_vector(2); 
    P_sat = @(Ai, Bi, Ci, bubbleT) antoine_equation(Ai, Bi, Ci, bubbleT);
    A_HCN = cmp(6).antaRT; % HCN
    B_HCN = cmp(6).antbRT; % HCN
    C_HCN = cmp(6).antcRT; % HCN
    A_H2O = cmp(1).antaRT; % H2O
    B_H2O = cmp(1).antbRT; % H2O
    C_H2O = cmp(1).antcRT; % H2O
    bT = x(1).*P_sat(A_H2O, B_H2O, C_H2O, bubbleT).*gamma_H2O+x(2).*P_sat(A_HCN, B_HCN, C_HCN, bubbleT).*gamma_HCN-P; % species1 = HCN, species2 = H2O when it comes to gamma; 
end

end