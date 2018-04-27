function [cmpout, untout, strout] = NH3_absorber_ideal(cmpin, untin, strin)
%%Main code for ammonia absorber

%% Defining variables

G = strin(5).G;
%G = 64;        %test value 
L = 250;
strout(7).G = G;
strout(7).L = L;

%molar fractions gaseous stream
y_HCN_in = strin(5).yHCN;
y_NH3_in = strin(5).yNH3;
y_H2_in = strin(5).yH2;
y_CH4_in = strin(5).yEgas;
y_N2_in = strin(5).yN2;
%y_HCN_in = 0.1706;
%y_NH3_in = 0.04;
%y_H2_in = 0.6876;
%y_CH4_in = 0.073;
%y_N2_in = 0.0288;


y_out = 1-(y_NH3_in-10^(-4));       %total gas stream out
y_HCN_out = y_HCN_in/y_out;
y_NH3_out = 10^(-4);                 %desired concentration
y_H2_out = y_H2_in/y_out;
y_CH4_out = y_CH4_in/y_out;
y_N2_out = y_N2_in/y_out;

strout(7).yHCN = y_HCN_out;
strout(7).yNH3 = y_NH3_out;
strout(7).yH2 = y_H2_out;
strout(7).yCH4 = y_CH4_out;
strout(7).yN2 = y_N2_out;

%molar fraction liquid stream
x_H2SO4_in = (G*y_NH3_in*1.1*0.5)/L;       
x_H2O_in = 1 - x_H2SO4_in;
%x_HCN_out = 0;                      %assumption, no physisorption
n_NH3 = G*(y_NH3_in-y_NH3_out);                 %assumption molar flow = mole fraction
n_H2SO4 = L*x_H2SO4_in;                %with or without a factor of 1.1?
n_tot = L*(x_H2O_in + x_H2SO4_in);
x_H2O_out = x_H2O_in;
x_H2SO4_out = (n_H2SO4 - 0.5*n_NH3)/n_tot;
x_ammoniumsulfate_out = x_H2SO4_in - x_H2SO4_out;

strout(6).xH2O = x_H2O_out;
strout(6).xH2SO4 = x_H2SO4_out;
strout(6).xAS = x_ammoniumsulfate_out;

%% Calculation for Heat Exchanger

%T_gas_reactor= 1600;
T_gas_reactor = strin(5).T;                 %temperature of inlet gas stream
T1_gas = 100+273;                    %desired value


deltaH_cooling_HCN = enthalpy_temperature(T_gas_reactor, cmpin, untin, 6) - enthalpy_temperature(T1_gas, cmpin, untin, 6);
deltaH_cooling_NH3 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 4) - enthalpy_temperature(T1_gas, cmpin, untin, 4);
deltaH_cooling_H2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 5) - enthalpy_temperature(T1_gas, cmpin, untin, 5);
deltaH_cooling_CH4 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 3) - enthalpy_temperature(T1_gas, cmpin, untin, 3);
deltaH_cooling_N2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 2) - enthalpy_temperature(T1_gas, cmpin, untin, 2);

Q_mixture = G * (deltaH_cooling_HCN*y_HCN_in + ...
                 deltaH_cooling_NH3*y_NH3_in +...
                 deltaH_cooling_H2*y_H2_in + ...
                 deltaH_cooling_CH4*y_CH4_in + ...
                 deltaH_cooling_N2*y_N2_in);

T_H2O_in = 15+273;
T_H2O_out = 90 + 273;             
deltaT1 = T1_gas - T_H2O_out;
deltaT2 = T1_gas - T_H2O_in;
             
LMTD = (deltaT1 - deltaT2)/(log(deltaT1/deltaT2));
             
A_exchanger = Q_mixture / (700*LMTD);
                 



%% First assumption: Physisorption neglected since chemisorption is mostly more effective, furthermore chemisorption is instantaneous and irreversible

% 2 NH3 + H2SO4 --> (NH4)2SO4               %reaction equation for absorption
% --> n_NH3 = 0.5 n_ammoniumsulfate         %moles


%calculation of temperature in ideal case, where heat capacity is not
%temperature dependent

%enthalphies for the incoming gas stream
hf_hcn = cmpin(6).deltaHf + cmpin(6).cp*(T1_gas-293);                   % Tobi fragen weil cmp(6).cp ergibt einen Vektor
hf_nh3 = cmpin(4).deltaHf + cmpin(4).cp*(T1_gas-293);
hf_h2 = cmpin(5).deltaHf + cmpin(5).cp*(T1_gas-293);
hf_ch4 = cmpin(3).deltaHf + cmpin(3).cp*(T1_gas-293);
hf_n2 = cmpin(2).deltaHf + cmpin(2).cp*(T1_gas-293);

%enthalpies for the outcoming gas stream
hf_hcn_out = @(T) cmpin(6).deltaHf + cmpin(6).cp*(T-293);           %not sure about the values in brackets
hf_nh3_out = @(T) cmpin(4).deltaHf + cmpin(4).cp*(T-293);
hf_h2_out = @(T) cmpin(5).deltaHf + cmpin(5).cp*(T-293);
hf_ch4_out = @(T) cmpin(3).deltaHf + cmpin(3).cp*(T-293);
hf_n2_out = @(T) cmpin(2).deltaHf + cmpin(2).cp*(T-293);

%enthalpies for the outcoming liquid stream
hf_h2so4_out = @(T) cmpin(7).deltaHf + cmpin(7).cp*(T-293);
hf_h2o_out = @(T) cmpin(1).deltaHf + cmpin(1).cp*(T-293);
hf_ammoniumsulfate = @(T) cmpin(8).deltaHf + cmpin(8).cp*(T-293);

%reaction enthalpy for the absorption of ammonia
h_rxn = 275020;                 %[J/mol]


E_gas_in = G.*(y_HCN_in.*(hf_hcn) + y_NH3_in.*(hf_nh3) + y_H2_in.*(hf_h2) + y_CH4_in.*(hf_ch4) + y_N2_in.*(hf_n2));
E_liquid_in = L.*(x_H2SO4_in*cmpin(7).deltaHf + x_H2O_in*cmpin(1).deltaHf);                 %no variation for enthalpies needed since the liquid has a temperature of 20 °C
E_gas_out = @(T) G.*(y_HCN_out*hf_hcn_out(T) + y_NH3_out*hf_nh3_out(T) + y_H2_out*hf_h2_out(T) + y_CH4_out*hf_ch4_out(T)+ y_N2_out*hf_n2_out(T));
E_liquid_out = @(T) L.*(x_H2SO4_out*hf_h2so4_out(T) + x_H2O_out*hf_h2o_out(T)+ x_ammoniumsulfate_out*hf_ammoniumsulfate(T));
E_rxn = G*h_rxn*(y_NH3_in - y_NH3_out)/2;

energy_balance = @(T) E_gas_in + E_liquid_in - E_gas_out(T) - E_liquid_out(T) + E_rxn;

T = fsolve(energy_balance, 300);
T_celsius = T-273;
%T=40+273;

strout(6).T = T;
strout(7).T = T;

%% Different approach


%% Calculation of the height of the column with HTU and NTU

%G = str(5).G;          % gas flow
velocity = 1;           %[m/s]
dia = diameter_calc2(velocity, G, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);                %randomly set diameter
r = dia/2;              %radius of column
A = r^2*pi;             %cross section of column
p = 101325;             %pressure in Pascal
%P_bar = 1;              % pressure in bar
a = 190;                %[m^2/m^3] specific surface of packing
ap = 492;               %[m^(-1)] Packing factor (dry)
dp = 0.025;             %[m] nominal diameter
R = 8.314472;           %universal gas constant [J/mol*K]




MW_NH3 = 17;            %molecular weights [g/mol]
MW_H2 = 2;
MW_HCN = 27;
MW_N2 = 28;
MW_CH4 = 16;
G_m = G*(y_NH3_in*MW_NH3 + y_H2_in*MW_H2 + y_HCN_in*MW_HCN + y_N2_in*MW_N2 + y_CH4_in*MW_CH4)/1000;         %[kg/s]

%assumption made that gaseous phase consists only of H2 and HCN
mu_g = mu_calc2(T);     %[kg/(m*s)]
rho = rho_calc(T, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in )/1000;      %[g/m^3]
%nu_g = mu_g/rho;
%Q = G*R*T/P_bar;

sum_V_H2 = 7.07;
sum_V_NH3 = 14.9;
D = 10^(-3)*T^1.75*(1/MW_H2 + 1/MW_NH3)^0.5/(((sum_V_H2^(1/3))+ (sum_V_NH3^(1/3)))^2)*10^(-4);          %[m^2/s]Frage mit 10^-4 als Faktor oder nicht?

Sc = mu_g/(rho*D);
%Sc = nu_g/D;


KG = 5.23*ap*D/(R*T) * (G_m/(ap*mu_g))^0.7 * Sc^(1/3) * (ap*dp)^(-2)*p;
%KG = 5.23*ap*D/R/T*(((Q/A)/(ap*nu_g))^0.7)*Sc^(1/3)*(ap*dp)^(-2)*p;
NTU = log(y_NH3_in/y_NH3_out);      %Colburn equation
HTU = G/(KG*a*A);
flow_ratio = L/G;

Z = NTU*HTU;
V_column = (dia/2)^2*pi*Z;

untout(2).htu = HTU;
untout(2).ntu = NTU;
untout(2).h = Z;
untout(2).V = V_column;


ratio = Z/dia;

[opex_H2O, opex_H2SO4, opex_wasterwater, opex_tot] = opex_calc(G, L,  x_H2O_in, x_H2SO4_in, x_H2O_out, x_H2SO4_out, x_ammoniumsulfate_out)
opex_heatexchanger = A_exchanger^0.65 * 25000 /10^6;        %[Mio. US$]
opex_tot2 = opex_tot+opex_heatexchanger;

capex = capex_calc(Z, dia)/10^6;

untout(2).capex = capex;
untout(2).opex = opex_tot2;

fprintf('Number of theoretical units: NTU = %g\n', NTU);
fprintf('Height of theoretical units: HTU = %g\n', HTU);
fprintf('Flow rate ratio: L/G = %g\n', flow_ratio);
fprintf('Outlet temperature [°C]: T = %g\n', T_celsius);
fprintf('Column Height [m]: H = %g\n', Z);
fprintf('Column Diameter [m]: D = %g\n', dia );
fprintf('NH3 Column CAPEX [Mio. US$]: Capex = %g\n', capex);
fprintf('NH3 Column OPEX [Mio. US$]: Opex = %g\n', opex_tot2);


end



