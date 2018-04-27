function [bubbleT] = bubblepoint(x, P, cmp, unt)
% calculates the bubble point for a binary mixture of H2O and HCN
% INPUT: x = vector of mole fractions; x(1) = H2O, x(2) = HCN; 
% P: pressure [Pa]
% cmp & unt: structs from Excel-file 



options = optimset('Display', 'off'); 
bubbleT_0 = 300; 
bubbleT = fsolve(@(bubbleT) bubbleT_solver(bubbleT), bubbleT_0, options); 


function [bT] = bubbleT_solver(bubbleT)
    delta_g12 = 1298.9610; % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1)
    delta_g21 = 539.9577; 
    alpha12 = 0.3836; % this has nothing to do with the relative volatility alpha 
    x1 = x(1); 
    gamma = nrtl(x1, bubbleT, delta_g12, delta_g21, alpha12); 
    gamma_HCN = gamma(1); 
    gamma_H2O = gamma(2); 
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