function [alpha] = relative_volatility(x1, T, cmpin, thermo_model)
% x1 = x_LK = x_HCN; 
A_H2O = 7.96681; B_H2O = 1668.21; C_H2O = 228.0; % (values for T > 60 °C), from Excel file  
A_HCN = 7.5282; B_HCN = 1329.5; C_HCN = 260.4;      % from Excel file
P_sat_H2O = @(temperature) antoine_equation_new(cmpin, temperature, 1); 
P_sat_HCN = @(temperature) antoine_equation_new(cmpin, temperature, 6);



delta_g12 = 1298.9610; % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1)
delta_g21 = 539.9577; 
alpha12 = 0.3836; % this has nothing to do with the relative volatility alpha

if strcmp(thermo_model, 'nrtl')
    gamma = @(TT) nrtl(x1, TT, delta_g12, delta_g21, alpha12); 
elseif strcmp(thermo_model, 'vanlaar')
    gamma = @(TT) vanlaar(x1, TT, delta_g12, delta_g21, alpha12); 
end
gamma_vector = gamma(T);
alpha = gamma_vector(1)*P_sat_HCN(T)/(gamma_vector(2)*P_sat_H2O(T));

  
    







end



