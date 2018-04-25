function [bubbleT] = bubblepoint(x, P, cmp, unt)
% calculates the bubble point for a binary mixture of H2O and HCN
% INPUT: x = vector of mole fractions; x(1) = H2O, x(2) = HCN; 
% P: pressure [Pa]
% cmp & unt: structs from Excel-file 



options = optimset('Display', 'off'); 
bubbleT = fsolve(@(bubbleT) bubbleT_solver, bubbleT_0, options); 


function [bT] = bubbleT_solver
    gamma = @(x1, T) nrtl(x1, T, delta_g12, delta_g21, alpha12); 
    gamma1 = gamma(1); 
    gamma2 = gamma(2); 
    P_sat = @(Ai, Bi, Ci, T) antoine_equation(Ai, Bi, Ci, T);
    A = cmp.antaRT; 
    B = cmp.antbRT; 
    C = cmp.antcRT;
    bT = x(1).*P_sat(A(1), B(1), C(1)).*gamma2+x(2).*P_sat(A(6), B(6), C(6)).*gamma2; % species1 = HCN, species2 = H2O when it comes to gamma; 
end

end