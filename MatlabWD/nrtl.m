function [gamma] = nrtl(x1, T, delta_g12, delta_g21, alpha12)
alpha21 = alpha12;                                                              % typical assumption in NRTL
R = 8.3144;                                                                     % ideal gas constant, [J/(mol*K)]
x2=1-x1;                                                                        % stoichiometric constraint
tau_12 = delta_g12/(R*T);                                                       % NRTL method, see https://en.wikipedia.org/wiki/Non-random_two-liquid_model 
tau_21 = delta_g21/(R*T); 
G12 = exp(-alpha12*tau_12); 
G21 = exp(-alpha21*tau_21); 
gamma(1) = exp(x2^2*(tau_21*((G21/(x1+x2*G21))^2)+tau_12*G12/((x2+x1*G12)^2)));   % activity coefficient of first species
gamma(2) = exp(x1^2*(tau_12*((G12/(x2+x1*G12))^2)+tau_21*G12/((x1+x2*G21)^2)));   % activity coefficient of second species
end



