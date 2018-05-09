function [cmpout, untout, strout] = NH3_absorber_ideal(cmpin, untin, strin)
%%Main code for ammonia absorber
cmpout=cmpin;
untout=untin;
strout=strin;
%% Defining variables

%Molar flows of gas and liquid inlet streams
G = strin(5).G;
L = 200;
strout(7).G = G;
strout(7).L = L;

%molar fractions inlet gas stream
y_HCN_in = strin(5).yHCN;
y_NH3_in = strin(5).yNH3;
y_H2_in = strin(5).yH2;
y_CH4_in = strin(5).yCH4;
y_N2_in = strin(5).yN2;

%Molar fractions outlet gas stream
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

%molar fractions inlet liquid stream
x_H2SO4_in = (G*y_NH3_in*1.1*0.5)/L;       
x_H2O_in = 1 - x_H2SO4_in;
%x_HCN_out = 0;                      %assumption, no physisorption
n_NH3 = G*(y_NH3_in-y_NH3_out);                 %assumption molar flow = mole fraction
n_H2SO4 = L*x_H2SO4_in;                %with or without a factor of 1.1?

%molar fractions outlet liquid stream
n_tot = L*(x_H2O_in + x_H2SO4_in);
x_H2O_out = x_H2O_in;
x_H2SO4_out = (n_H2SO4 - 0.5*n_NH3)/n_tot;
x_ammoniumsulfate_out = x_H2SO4_in - x_H2SO4_out;

strout(6).xH2O = x_H2O_out;
strout(6).xH2SO4 = x_H2SO4_out;
strout(6).xAS = x_ammoniumsulfate_out;

%Molar volumes for inlet stream
 V_m_h2so4 = 98.08 / 1.84;              %molecular weight divided by density [cm^3/mol]
 V_m_h2o = 18.01 / 0.9982;
 
%Volume flows of inlet stream
 V_h2so4 = x_H2SO4_in * L * V_m_h2so4;
 V_h2o = x_H2O_in * L * V_m_h2o;
 
%Volume percentage of sulfuric acid in the inlet liquid stream
 conc_H2SO4 = V_h2so4 / (V_h2so4 + V_h2o) *100;         


%% Calculation for Heat Exchanger

T_gas_reactor= 1600;
%T_gas_reactor = strin(5).T          %temperature of inlet gas stream
T1_gas = 100+273;                    %desired value

%Enthalpy difference of stream 5 for cooling down to 373K
deltaH_cooling_HCN = enthalpy_temperature(T_gas_reactor, cmpin, untin, 6) - enthalpy_temperature(T1_gas, cmpin, untin, 6);
deltaH_cooling_NH3 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 4) - enthalpy_temperature(T1_gas, cmpin, untin, 4);
deltaH_cooling_H2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 5) - enthalpy_temperature(T1_gas, cmpin, untin, 5);
deltaH_cooling_CH4 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 3) - enthalpy_temperature(T1_gas, cmpin, untin, 3);
deltaH_cooling_N2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 2) - enthalpy_temperature(T1_gas, cmpin, untin, 2);

%Heat flow of cooling down stream 5
Q_mixture = G * (deltaH_cooling_HCN*y_HCN_in + ...
                 deltaH_cooling_NH3*y_NH3_in +...
                 deltaH_cooling_H2*y_H2_in + ...
                 deltaH_cooling_CH4*y_CH4_in + ...
                 deltaH_cooling_N2*y_N2_in);

%Specifications for counter-current flowing cooling water
T_H2O_in = 15+273;
T_H2O_out = 90 + 273;             
deltaT1 = T1_gas - T_H2O_out;
deltaT2 = T1_gas - T_H2O_in;
             
LMTD = (deltaT1 - deltaT2)/(log(deltaT1/deltaT2));

%Area of heat exchanger, needed for Capex of heat exchanger
A_exchanger = abs(Q_mixture) / (700*LMTD);      %Q was negative, therefore a negative area was produced

%Mass flow of cooling water
mean_T = (T_H2O_in + T_H2O_out)/2;
cp_coolingwater = heat_capacity(mean_T, cmpin, untin, 1);
m_flow_water = Q_mixture * 18.01 /(cp_coolingwater * abs((T_H2O_in - T_H2O_out)))/1000/1000;        %in tons
                 



%% Assumption: Physisorption neglected since chemisorption is mostly more effective, furthermore chemisorption is instantaneous and irreversible

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
%h_rxn = 275020;                 %[J/mol]
h_rxn = cmpin(8).deltaHf - (cmpin(7).deltaHf + 2*cmpin(4).deltaHf);

%Energies of the streams calculated with the enthalpies
E_gas_in = G.*(y_HCN_in.*(hf_hcn) + y_NH3_in.*(hf_nh3) + y_H2_in.*(hf_h2) + y_CH4_in.*(hf_ch4) + y_N2_in.*(hf_n2));
E_liquid_in = L.*(x_H2SO4_in*cmpin(7).deltaHf + x_H2O_in*cmpin(1).deltaHf);                 %no variation for enthalpies needed since the liquid has a temperature of 20 �C
E_gas_out = @(T) G.*(y_HCN_out*hf_hcn_out(T) + y_NH3_out*hf_nh3_out(T) + y_H2_out*hf_h2_out(T) + y_CH4_out*hf_ch4_out(T)+ y_N2_out*hf_n2_out(T));
E_liquid_out = @(T) L.*(x_H2SO4_out*hf_h2so4_out(T) + x_H2O_out*hf_h2o_out(T)+ x_ammoniumsulfate_out*hf_ammoniumsulfate(T));
E_rxn = G*h_rxn*(y_NH3_in - y_NH3_out)/2;

%Calculation of outlet temperature with energy balance
energy_balance = @(T) E_gas_in + E_liquid_in - E_gas_out(T) - E_liquid_out(T)+ E_rxn;

options = optimset('Display','off');
T = fsolve(energy_balance, 300,options);
T_celsius = T-273;

strout(6).T = T;
strout(7).T = T;


%% Calculation of the height of the column with HTU and NTU

velocity = 1;           %[m/s]
dia = diameter_calc2(velocity, G, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);                
r = dia/2;              %radius of column
A = r^2*pi;             %cross section of column
p = 101325;             %pressure in Pascal
a = 190;                %[m^2/m^3] specific surface of packing
ap = 492;               %[m^(-1)] Packing factor (dry)
dp = 0.025;             %[m] nominal diameter
R = 8.314472;           %universal gas constant [J/mol*K]

%Molecular weights [g/mol]
MW_NH3 = 17;            
MW_H2 = 2;
MW_HCN = 27;
MW_N2 = 28;
MW_CH4 = 16;
G_m = G*(y_NH3_in*MW_NH3 + y_H2_in*MW_H2 + y_HCN_in*MW_HCN + y_N2_in*MW_N2 + y_CH4_in*MW_CH4)/1000;         %[kg/s]


%Calculation for Diffusivity
sum_V_H2 = 7.07;
sum_V_NH3 = 14.9;
D = 10^(-3)*T^1.75*(1/MW_H2 + 1/MW_NH3)^0.5/(((sum_V_H2^(1/3))+ (sum_V_NH3^(1/3)))^2)*10^(-4);          %[m^2/s]Frage mit 10^-4 als Faktor oder nicht?

%Calculation of Schmidt number
%assumption made that gaseous phase consists only of H2 and HCN
mu_g = mu_calc2(T);                                                             %[kg/(m*s)]
rho = rho_calc(T, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in )/1000;        %[g/m^3]
Sc = mu_g/(rho*D);

%Calculation of mass transfer coefficient for gas 
KG = 5.23*ap*D/(R*T) * (G_m/(ap*mu_g))^0.7 * Sc^(1/3) * (ap*dp)^(-2)*p;

NTU = log(y_NH3_in/y_NH3_out);      %Colburn equation
HTU = G/(KG*a*A);



%Calculation of height und volume of column with consideration of free
%volume share in the column
Z = NTU*HTU;
free_vol = 0.74;                            %free volume share in the column
V_column = (dia/2)^2*pi*Z / free_vol;
dia_true = (V_column /(pi*Z))^0.5 *2;       %diameter with the consideration of free volume share in column

untout(2).dia = dia_true;

%Calculating ratios, ratio of height to diameter should be between 5 and 15
flow_ratio = L/G;
ratio = Z/dia_true;

untout(2).ratio = flow_ratio;

%Calculation of Opex
[opex_H2O, opex_H2SO4, opex_wasterwater, opex_tot] = opex_calc(G, L,  x_H2O_in, x_H2SO4_in, x_H2O_out, x_H2SO4_out, x_ammoniumsulfate_out);     
time = 8000*60*60;
Opex_cooling_water = m_flow_water*time*0.10;              %[US$]

opex_tot2 = opex_tot + Opex_cooling_water;
    %assigning values to get rid of warnings
     opex_H2O=opex_H2O;
     opex_H2SO4 = opex_H2SO4;
     opex_wastewater = opex_wasterwater;

%Calculation of Capex, including heat exchanger
capex = capex_calc(Z, dia);
capex_heatexchanger = A_exchanger^0.65 * 25000;  
capex_tot = capex_heatexchanger + capex;


%VARS OUT

untout(2).htu = HTU;    %height theoretical units NH3 ADSORBER
untout(2).ntu = NTU;    %number theoritac units
untout(2).h = Z;        %height NH3 Absorber
untout(2).V = V_column; %VOlume
untout(2).capex = capex_tot;
untout(2).opex = opex_tot2;

% fprintf('Number of theoretical units: NTU = %g\n', NTU);
% fprintf('Height of theoretical units: HTU = %g\n', HTU);
% fprintf('Flow rate ratio: L/G = %g\n', flow_ratio);
% fprintf('Outlet temperature [�C]: T = %g\n', T_celsius);
% fprintf('Column Height [m]: H = %g\n', Z);
% fprintf('Column Diameter [m]: D = %g\n', dia );
% fprintf('NH3 Column CAPEX [Mio. US$]: Capex = %g\n', capex);
% fprintf('NH3 Column OPEX [Mio. US$]: Opex = %g\n', opex_tot2);


end



