function [cmpout, untout, strout] = NH3_absorber_nonideal4(cmpin, untin, strin)
%% Modelling for Ammonia Absorber with a temperature dependent heat capacity
cmpout=cmpin;
untout=untin;
strout=strin;
%% Defining variables

G = strin(5).G;
%G = 64;                             %test value 
L = 100;                            %molar fraction liquid stream
strout(7).G = G;
strout(7).L = L;

%molar fractions inlet gaseous stream
y_HCN_in = strin(5).yHCN;
y_NH3_in = strin(5).yNH3;
y_H2_in = strin(5).yH2;
y_CH4_in = strin(5).yEgas;
y_N2_in = strin(5).yN2;

%y_HCN_in = 0.205;
%y_NH3_in = 0.0128;
%y_H2_in = 0.6876;
%y_CH4_in = 0.073;
%y_N2_in = 0.0288;

%molar fractions of outlet gas stream
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

%molar fractions of inlet liquid stream
x_H2SO4_in = (G*y_NH3_in*1.2*0.5)/L;       
x_H2O_in = 1 - x_H2SO4_in;
%y_HCN_out = 0;                                 %assumption, no physisorption
n_NH3 = G*(y_NH3_in-y_NH3_out);                 %assumption molar flow = mole fraction
n_H2SO4 = L*x_H2SO4_in;                

%molar fractions of outlet liquid stream
n_tot = L*(x_H2O_in + x_H2SO4_in);
x_H2O_out = x_H2O_in;
x_H2SO4_out = (n_H2SO4 - 0.5*n_NH3)/n_tot;
x_ammoniumsulfate_out = (x_H2SO4_in - x_H2SO4_out);

strout(6).xH2O = x_H2O_out;
strout(6).xH2SO4 = x_H2SO4_out;
strout(6).xAS = x_ammoniumsulfate_out;

%molar volumes
 V_m_h2so4 = 98.08/1.84;            %molecular weights divided by the densities                   
 V_m_h2o = 18.01/0.9982;
 
 %Volume of liquid inlet components
 V_h2so4 = x_H2SO4_in*L*V_m_h2so4;              
 V_h2o = x_H2O_in * L * V_m_h2o;
 
 %concentration of sulfuric acid in volume percent
 conc_H2SO4 = V_h2so4/(V_h2so4 + V_h2o)*100;            


%% Calculation for Heat Exchanger

T_gas_reactor= 1600;
%T_gas_reactor = strin(5).T;                    %temperature of inlet gas stream
T1_gas = 100+273;                               %desired value

%Enthalpies of components of stream 5 for cooling
deltaH_cooling_HCN = enthalpy_temperature(T_gas_reactor, cmpin, untin, 6) - enthalpy_temperature(T1_gas, cmpin, untin, 6);
deltaH_cooling_NH3 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 4) - enthalpy_temperature(T1_gas, cmpin, untin, 4);
deltaH_cooling_H2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 5) - enthalpy_temperature(T1_gas, cmpin, untin, 5);
deltaH_cooling_CH4 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 3) - enthalpy_temperature(T1_gas, cmpin, untin, 3);
deltaH_cooling_N2 = enthalpy_temperature(T_gas_reactor, cmpin, untin, 2) - enthalpy_temperature(T1_gas, cmpin, untin, 2);

%Heat flow
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
    
%Area of heat exchanger needed for Capex
A_exchanger = abs(Q_mixture) / (700*LMTD);

%Mass flow of cooling water
mean_T = (T_H2O_in + T_H2O_out)/2;
cp_coolingwater = heat_capacity(mean_T, cmpin, untin, 1);
m_flow_water = Q_mixture * 18.01 /(cp_coolingwater * (T_gas_reactor - T_H2O_out))/1000/1000;        %in tons


%% Assumption: Physisorption neglected since chemisorption is mostly more effective, furthermore chemisorption is instantaneous and irreversible

% 2 NH3 + H2SO4 --> (NH4)2SO4               %reaction equation for absorption
% --> n_NH3 = 0.5 n_ammoniumsulfate         %moles

%specification of temperatures and reducing them for shomate equation
T0_liquid = 20+273;                         %input temperature of liquid in [K]
T0_gas = 100+273;
t0_liquid = T0_liquid/1000;                 %reduced temperature for Shomate equation
t0_gas = T0_gas/1000;


% Shomate coefficients
%[HCN, NH3, H2, CH4, N2, H2SO4, H2O, ammoniumsulfate]
A = [32.69, 19.99, 33.06, -0.703, 19.50, 47.28, -203.60];
B = [22.59, 49.77, -11.36, 108.47, 19.88, 190.33, 1523.29];
C = [-4.37, -15.38, 11.43, -42.52, -8.59, -148.13, -3196.41];
D = [-0.407, 1.92, -2.77, 5.86, 1.36, 43.86, 2474.45];
E = [-0.28, 0.19, -0.16, 0.678, 0.527, -0.74, 3.855];
F = [123.5, -53.31, -9.98, -76.84, -4.93, -758.95, -256.54];
H = [135.14, -45.89, 0.0, -74.87, 0.0, -735.12, -285.82];

% Enthalpies for the incoming gas stream
hf_hcn_non = cmpin(6).deltaHf + ((A(1)*t0_gas) + (B(1)*((t0_gas)^2)/2) + (C(1)*((t0_gas)^3)/3) + (D(1)*((t0_gas)^4)/4) - (E(1)/t0_gas) +F(1) - H(1))*1000;                   
hf_nh3_non = cmpin(4).deltaHf + ((A(2)*t0_gas) + (B(2)*((t0_gas)^2)/2) + (C(2)*((t0_gas)^3)/3) + (D(2)*((t0_gas)^4)/4) - (E(2)/t0_gas) +F(2) - H(2))*1000;
hf_h2_non = cmpin(5).deltaHf + ((A(3)*t0_gas) + (B(3)*((t0_gas)^2)/2) + (C(3)*((t0_gas)^3)/3) + (D(3)*((t0_gas)^4)/4) - (E(3)/t0_gas) +F(3) - H(3))*1000; 
hf_ch4_non = cmpin(3).deltaHf + ((A(4)*t0_gas) + (B(4)*((t0_gas)^2)/2) + (C(4)*((t0_gas)^3)/3) + (D(4)*((t0_gas)^4)/4) - (E(4)/t0_gas) +F(4) - H(4))*1000; 
hf_n2_non = cmpin(2).deltaHf + ((A(5)*t0_gas) + (B(5)*((t0_gas)^2)/2) + (C(5)*((t0_gas)^3)/3) + (D(5)*((t0_gas)^4)/4) - (E(5)/t0_gas) +F(5) - H(5))*1000; 

% Enthalpie for incoming liquid stream
hf_h2o_non = cmpin(1).deltaHf + (((A(7)*t0_liquid) + (B(7)*((t0_liquid)^2)/2) + (C(7)*((t0_liquid)^3)/3) + (D(7)*((t0_liquid)^4)/4) - (E(7)/t0_liquid) +F(7) - H(7))*1000);
hf_h2so4_non =  cmpin(7).deltaHf + (((A(6)*t0_liquid) + (B(6)*((t0_liquid)^2)/2) + (C(6)*((t0_liquid)^3)/3) + (D(6)*((t0_liquid)^4)/4) - (E(6)/t0_liquid) +F(6) - H(6))*1000);

%enthalpies for the outcoming gas stream
hf_hcn_out_non = @(t_out) cmpin(6).deltaHf +(((A(1)*t_out) + (B(1)*((t_out)^2)/2) + (C(1)*((t_out)^3)/3) + (D(1)*((t_out)^4)/4) - (E(1)/t_out) +F(1) - H(1))*1000);  
hf_nh3_out_non = @(t_out) cmpin(4).deltaHf + (((A(2)*t_out) + (B(2)*((t_out)^2)/2) + (C(2)*((t_out)^3)/3) + (D(2)*((t_out)^4)/4) - (E(2)/t_out) +F(2) - H(2))*1000); 
hf_h2_out_non = @(t_out) cmpin(5).deltaHf + (((A(3)*t_out) + (B(3)*((t_out)^2)/2) + (C(3)*((t_out)^3)/3) + (D(3)*((t_out)^4)/4) - (E(3)/t_out) +F(3) - H(3))*1000);  
hf_ch4_out_non = @(t_out) cmpin(3).deltaHf + (((A(4)*t_out) + (B(4)*((t_out)^2)/2) + (C(4)*((t_out)^3)/3) + (D(4)*((t_out)^4)/4) - (E(4)/t_out) +F(4) - H(4))*1000);  
hf_n2_out_non = @(t_out) cmpin(2).deltaHf + (((A(5)*t_out) + (B(5)*((t_out)^2)/2) + (C(5)*((t_out)^3)/3) + (D(5)*((t_out)^4)/4) - (E(5)/t_out) +F(5) - H(5))*1000);  

%enthalpies for the outcoming liquid stream
hf_h2o_out_non = @(t_out) cmpin(1).deltaHf + (((A(7)*t_out) + (B(7)*((t_out)^2)/2) + (C(7)*((t_out)^3)/3) + (D(7)*((t_out)^4)/4) - (E(7)/t_out) +F(7) - H(7))*1000);  
hf_h2so4_out_non = @(t_out) cmpin(7).deltaHf + (((A(6)*t_out) + (B(6)*((t_out)^2)/2) + (C(6)*((t_out)^3)/3) + (D(6)*((t_out)^4)/4) - (E(6)/t_out) +F(6) - H(6))*1000);  
hf_ammoniumsulfate_non = @(t_out) cmpin(8).deltaHf + cmpin(8).cp*(t_out-293/1000);

%reaction enthalpy for the absorption of ammonia
h_rxn = 275020;                 %[J/mol]

%Energies of the streams and of the reaction with the enthalpies
E_gas_in_non = G.*(y_HCN_in.*(hf_hcn_non) + y_NH3_in.*(hf_nh3_non) + y_H2_in.*(hf_h2_non) + y_CH4_in.*(hf_ch4_non) + y_N2_in.*(hf_n2_non));
E_liquid_in_non = L.*(x_H2SO4_in*hf_h2so4_non + x_H2O_in*hf_h2o_non);                 
E_gas_out_non = @(t_out) G.*(y_HCN_out*hf_hcn_out_non(t_out) + y_NH3_out*hf_nh3_out_non(t_out) + y_H2_out*hf_h2_out_non(t_out) + y_CH4_out*hf_ch4_out_non(t_out)+ y_N2_out*hf_n2_out_non(t_out));
E_liquid_out_non = @(t_out) L.*(x_H2SO4_out*hf_h2so4_out_non(t_out) + x_H2O_out*hf_h2o_out_non(t_out)+ x_ammoniumsulfate_out*hf_ammoniumsulfate_non(t_out));
E_rxn = G*h_rxn*(y_NH3_in - y_NH3_out)/2;

energy_balance_non = @(t_out) (E_gas_in_non + E_liquid_in_non - E_gas_out_non(t_out) - E_liquid_out_non(t_out) + E_rxn);

%t_non = lsqnonlin(energy_balance_non, 50);
t_non = fsolve(energy_balance_non, 0.3);
T_non = t_non*1000;                  %[K]

strout(6).T = T_non;
strout(7).T = T_non;

%% Calculation of the height of the column with HTU and NTU

velocity = 1;           %[m/s]
dia = diameter_calc2(velocity, G, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);                %randomly set diameter
r = dia/2;              %radius of column
A = r^2*pi;             %cross section of column
p = 101325;             %pressure in Pascal
a = 190;                %[m^2/m^3] specific surface of packing
ap = 492;               %[m^(-1)] Packing factor (dry)
dp = 0.025;             %[m] nominal diameter
R = 8.314472;           %universal gas constant [J/mol*K]



%molecular weights [g/mol]
MW_NH3 = 17;            
MW_H2 = 2;
MW_HCN = 27;
MW_N2 = 28;
MW_CH4 = 16;
G_m = G*(y_NH3_in*MW_NH3 + y_H2_in*MW_H2 + y_HCN_in*MW_HCN + y_N2_in*MW_N2 + y_CH4_in*MW_CH4)/1000;         %[kg/s]




%Calculation of Diffusivity
sum_V_H2 = 7.07;
sum_V_NH3 = 14.9;
D = 10^(-3)*T_non^1.75*(1/MW_H2 + 1/MW_NH3)^0.5/(((sum_V_H2^(1/3))+ (sum_V_NH3^(1/3)))^2)*10^(-4);          %[m^2/s]Frage mit 10^-4 als Faktor oder nicht?

%Calculation of Schmidt number
%assumption made that gaseous phase consists only of H2 and HCN
mu_g = mu_calc2(T_non);     %[kg/(m*s)]
rho = rho_calc(T_non, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in )/1000;      %[g/m^3]
Sc = mu_g/(rho*D);

%Calculation of mass transfer coefficient for gas
KG = 5.23*ap*D/(R*T_non) * (G_m/(ap*mu_g))^0.7 * Sc^(1/3) * (ap*dp)^(-2)*p;

NTU = log(y_NH3_in/y_NH3_out);      %Colburn equation
HTU = G/(KG*a*A);

%Calculation of height and volume of column, considering the free volume
%share of the volume
Z = NTU*HTU;
free_vol = 0.74;                            %free volume share in the column
V_column = (dia/2)^2*pi*Z / free_vol;
dia_true = (V_column /(pi*Z))^0.5 *2;       %diameter with the consideration of free volume share in column

flow_ratio = L/G;
ratio = Z/dia_true;

%Calculation for Opex
[opex_H2O, opex_H2SO4, opex_wasterwater, opex_tot] = opex_calc(G, L,  x_H2O_in, x_H2SO4_in, x_H2O_out, x_H2SO4_out, x_ammoniumsulfate_out);
time = 8000*60*60;
Opex_cooling_water = m_flow_water*time*0.10;              %[US$]
opex_tot2 = opex_tot + Opex_cooling_water;
    %assigning values to get rid of warnings
     opex_H2O=opex_H2O;
     opex_H2SO4 = opex_H2SO4;
     opex_wastewater = opex_wasterwater;

%Calculation for Capex including heat exchanger
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
% fprintf('Outlet temperature [K]: T = %g\n', T_non);
% fprintf('Column Height [m]: H = %g\n', Z);
% fprintf('Column Diameter [m]: D = %g\n', dia );
% fprintf('NH3 Column CAPEX [Mio. US$]: Capex = %g\n', capex_tot);
% fprintf('NH3 Column OPEX [Mio. US$]: Opex = %g\n', opex_tot2);

 
 end  


