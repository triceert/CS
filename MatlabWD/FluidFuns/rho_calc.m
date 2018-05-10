function rho = rho_calc(t, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in );

%assumption of ideal gas for density calculation

T = [293, 303, 313, 323, 333, 343, 353];
rho_general = [0.041, 0.039, 0.038, 0.037, 0.036, 0.035, 0.034].*1000;        %[mol/m^3] calculated from the ideal gas law

rho_HCN = rho_general*27;               %[g/m^3]
rho_H2 = rho_general*2;
rho_NH3 = rho_general*17;
rho_CH4 = rho_general*16;
rho_N2 = rho_general*28;

rho_mix = rho_HCN*y_HCN_in + rho_H2*y_H2_in + rho_NH3*y_NH3_in + rho_CH4*y_CH4_in + rho_N2*y_N2_in;            %combined densities but without temperature dependence

rho_mix_temp = fitlm(T, rho_mix);

intercept = table2array(rho_mix_temp.Coefficients(1,1));
slope = table2array(rho_mix_temp.Coefficients(2,1));

rho = intercept + slope*t;

end