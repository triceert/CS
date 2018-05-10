%  HCN ABSORPTION IDEAL CALCULATION

function [cmpout, untout, strout] = hcnideal(cmpin, untin, strin)
cmpout=cmpin;
untout=untin;
strout=strin;

% General assumption: Isothermal temperature profile

% Gas flow inlet

G_in = strin(7).G;             % Inlet gas molar stream [mol/s]  
 
y_HCN_in = strin(7).yHCN;      % Molar fractions in the gas stream
y_H2_in = strin(7).yH2;
y_CH4_in = strin(7).yCH4;
y_NH3_in = strin(7).yNH3;
y_N2_in = strin(7).yN2;

y_in = [y_HCN_in, y_NH3_in, y_CH4_in,  y_H2_in,  y_N2_in];  % molar fractions of the gas in an array

G_HCN_in = G_in * y_HCN_in;    % Molar flows of each species [mol/s]
G_H2_in = G_in * y_H2_in;
G_CH4_in = G_in * y_CH4_in;
G_NH3_in = G_in * y_NH3_in;
G_N2_in = G_in * y_N2_in;

G_in_i = [G_NH3_in, G_HCN_in, G_CH4_in, G_H2_in, G_N2_in]; % Input Molar Flows as array [NH3 HCN CH4 H2 N2]

% Gas flow outlet

G_out= -sum(G_in_i(1)+G_in_i(3:end))/(10^(-4)-1);          % Total output molar flow [mol/s]
strout(8).G = G_out;                                       % Passes the output molar flow to the Excel-file

y_HCN_out = 1e-4;                                          % Outlet molar fractions of the diff species
y_NH3_out = G_in_i(1)/G_out;   
y_CH4_out = G_in_i(3)/G_out; 
y_H2_out  = G_in_i(4)/G_out;
y_N2_out  = G_in_i(5)/G_out;
y_out = [y_HCN_out, y_NH3_out, y_CH4_out, y_H2_out, y_N2_out];

strout(8).yHCN = y_HCN_out;
strout(8).yNH3 = y_NH3_out;
strout(8).yCH4 = y_CH4_out;
strout(8).yH2 = y_H2_out;
strout(8).yN2 = y_N2_out;

% Inlet Liquid flow

L_in = strin(10).L;       % Input liquid total molar stream from distillation, calculated as 1.5*L(min) [mol/s]
L_out= L_in+(G_in_i(2)-(G_out*10^(-4)));     % Output total liquid molar stream in [mol/s]
strout(9).L = L_out;   

% Made assumption: No HCN in the inlet liquid stream 

x_HCN_in = 0;
x_H2O_in = 1;
x_in = [x_HCN_in, x_H2O_in];   

% Outlet liquid flow

x_HCN_out = (G_in_i(2) - y_HCN_out*G_out)/(L_out);    % Outlet liquid molar fraction
x_H2O_out = 1-x_HCN_out;
x_out = [x_HCN_out, x_H2O_out];

strout(9).xHCN = x_HCN_out;
strout(9).xH2O = x_H2O_out;

G_in_i = G_in.*y_in;                   % Inlet and outlet molar streams of each species
G_out_i= G_out.*y_out;
L_in_i = L_in.*x_in;
L_out_i= L_out.*x_out;

% Inlet temperature of the gas and liquid flow
T_gas_in = strin(7).T; 
T_liquid_in = 273 + 20;              % Temperature of the inlet liquid as said in the assignement
t_gas_in = T_gas_in/1000;            % Modified temperature used for the Shomate equations
t_liquid_in = T_liquid_in/1000;

% Standard formation enthalpy [J/mol] (from literature)
h0f_HCN = cmpin(6).deltaHf ;
h0f_H2O = cmpin(1).deltaHf ;
h0f_H2 =  cmpin(5).deltaHf ;
h0f_N2 =  cmpin(2).deltaHf ;
h0f_NH3 = cmpin(4).deltaHf ;
h0f_CH4 = cmpin(3).deltaHf ;

% Heat capacities at 298.15 K and standard pressure [J/mol K] (from literature)
cp_0_HCN = cmpin(6).cp ;      
cp_0_H2O = cmpin(1).cp ;
cp_0_H2 =  cmpin(5).cp ;
cp_0_N2 = cmpin(2).cp;
cp_0_NH3 = cmpin(4).cp;
cp_0_CH4 = cmpin(3).cp ;

% Made assumption 2: Heat capacitites of the species are not temperature dependant
% Modified formation enthalpies using hf = h0f + cp (T1-T0)

% Enthalpies of the inlet gas stream [J/mol]
hf_HCN_gas_in =  h0f_HCN + (cp_0_HCN * (T_gas_in - 293.15));  
hf_H2_gas_in =  h0f_H2 + (cp_0_H2 * (T_gas_in -  293.15)); 
hf_N2_gas_in =  h0f_N2 + (cp_0_N2 *(T_gas_in -  293.15)); 
hf_NH3_gas_in =  h0f_NH3 + (cp_0_NH3*  (T_gas_in -  293.15)); 
hf_CH4_gas_in = h0f_CH4 + (cp_0_CH4 *(T_gas_in -  293.15)); 

% Enthalpies of the inlet liquid stream [J/mol]

hf_H2O_liquid_in = h0f_H2O + (cp_0_H2O *(T_liquid_in -  293.15));

% Enthalpies of the outlet gas stream [J/mol]

hf_HCN_gas_out = @(T_out) h0f_HCN + (cp_0_HCN * (T_out -  293.15));  
hf_H2_gas_out = @(T_out) h0f_H2 + (cp_0_H2 * (T_out -  293.15)); 
hf_N2_gas_out = @(T_out) h0f_N2 + (cp_0_N2 * (T_out - 293.15)); 
hf_NH3_gas_out = @(T_out) h0f_NH3 + (cp_0_NH3 * (T_out-  293.15)); 
hf_CH4_gas_out = @(T_out) h0f_CH4 + (cp_0_CH4 * (T_out-  293.15));

% Enthalpies of the outlet liquid stream [J/mol]

hf_HCN_liquid_out = @(T_out) h0f_HCN + (cp_0_HCN * (T_out-  293.15));  
hf_H20_liquid_out = @(T_out) h0f_H2O + (cp_0_H2O * (T_out -  293.15)); 

% Assumption 3: Heat of absorption not considered for the energy balance
% Energy balance 

E_gas_in = G_in * ((y_HCN_in * hf_HCN_gas_in) + ...
            (y_H2_in * hf_H2_gas_in) + ...
            (y_N2_in * hf_N2_gas_in) + ...
            (y_NH3_in * hf_NH3_gas_in) + ...
            (y_CH4_in * hf_CH4_gas_in)); 
        
E_liquid_in = L_in * (x_H2O_in * hf_H2O_liquid_in);

E_gas_out = @(T_out) G_out * ((y_HCN_out * hf_HCN_gas_out(T_out)) + ...
            (y_H2_out * hf_H2_gas_out(T_out)) + ...
            (y_N2_out * hf_N2_gas_out(T_out)) + ...
            (y_NH3_out * hf_NH3_gas_out(T_out)) + ...
            (y_CH4_out * hf_CH4_gas_out(T_out)));
        
E_liquid_out = @(T_out) L_out * ((x_H2O_out * hf_H20_liquid_out(T_out)) + ...
            (x_HCN_out * hf_HCN_liquid_out(T_out))); 
 
Ebalance = @(T_out) E_gas_in + E_liquid_in - E_gas_out(T_out) - E_liquid_out(T_out);
optis=optimset('Display','off');
Temp_out = fsolve (Ebalance, 300, optis);     % Outlet temperature [K]
Temp_out_celsius = Temp_out - 273;     % Outlet temperature [�C]
strout(9).T = Temp_out;
strout(8).T = Temp_out;

flow_ratio = L_in/G_in;         
flow_ratio_out = L_out/G_out;

untout(3).ratio = flow_ratio;

 %%
 % Calculation of the height of the HCN absorber column (ideal)
 
M_H2O = 18;    % Molecular weight of the different species [g/mol]
M_HCN = 27;
M_H2 = 2;
M_N2 = 28;
M_NH3 = 17;
M_CH4 = 16;
 
% Given values

p = 101325;      % Pressure inside the column [Pa]
P_bar = 1.013;   % Pressure [Bar]
p_atm = 1;       % Pressure [atm]
velocity = 2;    % Gas velocity in [m/s]
dia = diameter_calc2 (velocity,G_in, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);    % Diameter in m
area = pi*((dia/2)^2);       % Cross sectional area of the column [m^2] 
v_frei = 0.74;               % Free volume fraction
tot_surf = 190;              % Total surface in m^2/m^3
ap = 492;                    % Packing factor (dry) [m^(-1)]
aw = 587;                    % Packing factor (wet) [m^(-1)]
g = 9.81;                    % Acceleration of gravity [m/s^2]
dp = 0.025;                  % Nominal diameter [m]
R = 8.314;                   % Universal gas constant [J/mol*K]


% Assumption 4: Resistivity in the liquid phase neglected
% Assumption 5: Gas phase consists only of H2 and HCN

% HG Value

mu_sum = mu_calc2(Temp_out);     % Dynamic viscosity [kg/(m*s)]
G_m = G_in*(y_NH3_in*M_NH3 + y_H2_in*M_H2 + y_HCN_in*M_HCN + y_N2_in*M_N2 + y_CH4_in*M_CH4)/1000;  %[kg/s]
rho_sum = rho_calc(Temp_out, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in)/1000;   % Densitites in dependance of the outlet temperature [g/m^3]

sum_V_H2 = 7.07;      % from literature
sum_V_HCN = 24.17  ;
D = 10^(-3)*Temp_out^1.75*(1/M_H2 + 1/M_HCN)^0.5/((((sum_V_H2^(1/3))+ (sum_V_HCN^(1/3)))^2)*p_atm)*10^(-4);    % Diffusion coefficient   
Sc = mu_sum/(rho_sum*D);   % Schmidt number
KG = 5.23*ap*D/(R*Temp_out) * (G_m/(ap*mu_sum))^0.7 * Sc^(1/3) * (ap*dp)^(-2)*p;
HG = G_in/(area * KG * tot_surf);
HTU = HG;            % Height of theoretical unit [m]                                  

henry_HCN = HenrysConstant(Temp_out,cmpin,6) ;    % Henry coefficient of HCN in [Pa]
m = henry_HCN/p;                                  % Equilibrium constant
A = L_in/(G_in * m);
alpha = (y_HCN_in-y_HCN_out)/(y_HCN_in -m*x_HCN_in);
NTU= A./(A-1).*log((1-alpha./A)/(1-alpha));       % Number of theoretical units from SPt
h = HTU * NTU;                                    % Height of the absorber [m]
ratio = h/dia;                                     % Should be ideally between 5 and 15
V_column = (dia/2)^2 * h * pi/0.74;                    % Volume of the column [m^3]
dia_true = 2 * (V_column/(h * pi))^0.5;           % True diameter with considering the availalbe area

untout(3).dia = dia_true;
%%
% CAPEX calculation

CAPEX_hcnabs = capex_calc(h, dia_true);
CAPEX_mil = CAPEX_hcnabs/1000000;

% OPEX calculation

[opex_H2O, opex_H2SO4, opex_wastewater, opex_tot] = opex_calc(G_in, L_in, x_H2O_in, 0, 0, 0, 0);


% Calculation of the minimum liquid inlet stream required

min_ratio = ((y_HCN_in - y_HCN_out)/((y_HCN_in/m)-x_HCN_in));
L_min = G_in * min_ratio;
L_true = L_min * 1.5;         % Real liquid stream , factor 1.5 randomly picked            

%fprintf('Number of theoretical units: NTU = %g\n', NTU);
%fprintf('Height of theoretical units: HTU = %g\n', HTU);
%fprintf('Flow rate ratio: L/G = %g\n', flow_ratio);
%fprintf('Outlet temperature [�C]: T = %g\n', Temp_out_celsius);
%fprintf('Column Height [m]: H = %g\n', h);
%fprintf('Column Diameter [m]: D = %g\n', dia );
%fprintf('HCN Column CAPEX [Mio. US$]: Capex = %g\n', CAPEX_mil);
%fprintf('HCN Column OPEX [Mio. US$]: Opex = %g\n', opex_tot);

untout(3).h = h;
untout(3).capex = CAPEX_hcnabs;
untout(3).opex = opex_tot;
untout(3).ntu = NTU;
untout(3).htu = HTU;
untout(3).V = V_column;

strout(9).p = p;   % pressure in Pa 
end












 
 
 
 
 
 



 
 


