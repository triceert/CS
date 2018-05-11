function mu_g = mu_calc2(t)

%Assumptions made: Binary gas phase consisting of 20% hcn and 80% hydrogen
%                  linear correlation of mu und temperature

T = [293, 303, 313, 323, 333, 343, 353];
mu_H2 = [8.745, 8.952, 9.154, 9.353, 9.549, 9.742, 9.929].*10^(-6);    %Source: https://www.lmnoeng.com/Flow/GasViscosity.php
mu_HCN = [17.419, 17.864, 18.302, 18.734, 19.160, 19.580, 19.994].*10^(-6);     %Source: matlab code
mu_mix = 0.2.*mu_HCN + 0.8*mu_H2;

lin_mod = fitlm(T, mu_mix);

intercept = table2array(lin_mod.Coefficients(1,1));
slope = table2array(lin_mod.Coefficients(2,1));

mu_g = intercept + slope*t;



end
