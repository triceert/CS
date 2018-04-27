% HCN Absorption column all

function [cmpout, untout, strout] = hcnabsorption2(cmpin, untin, strin)


%G_in = 64;             % Molar gas stream [mol/s]
 G_in = strin(7).G;
 
 y_HCN_in = strin(7).yHCN;      % Molar fractions in the gas stream
 y_H2_in = strin(7).yH2;
 y_CH4_in = strin(7).yCH4;
 y_NH3_in = strin(7).yNH3;
 y_N2_in = strin(7).yN2;
 
 %y_HCN_in = 0.2076;
 %y_H2_in = 0.6964;
 %y_CH4_in = 0.0739;
 %y_NH3_in = 0.0001;
 %y_N2_in = 0.0292;
 
 G_HCN_in = G_in * y_HCN_in;           % Molar flows of each species
 G_H2_in = G_in * y_H2_in;
 G_CH4_in = G_in * y_CH4_in;
 G_NH3_in = G_in * y_NH3_in;
 G_N2_in = G_in * y_N2_in;
 
 G_in_i = [G_NH3_in, G_HCN_in, G_CH4_in, G_H2_in, G_N2_in]; %Input Molar Flows as array [NH3 HCN CH4 H2 N2]
 G_out= -sum(G_in_i(1)+G_in_i(3:end))/(10^(-4)-1);          % Output total molar flow
 strout(8).G = G_out;
 
 L_in = 518.1;                                % Input liquid total molar stream from distillation
 L_out= L_in+(G_in_i(2)-(G_out*10^(-4)));     % Output total liquid molar stream in [mol/s]
 strout(9).L = L_out;
y_in = [y_HCN_in, y_NH3_in, y_CH4_in,  y_H2_in,  y_N2_in];  %molal fractions of the gas in an array

y_HCN_out = 1e-4;                 
y_NH3_out = 1e-4;   
y_CH4_out = G_in_i(3)/G_out; 
y_H2_out  = G_in_i(4)/G_out;
y_N2_out  = G_in_i(5)/G_out;
y_out = [y_HCN_out, y_NH3_out, y_CH4_out, y_H2_out, y_N2_out];

 strout(8).yHC = y_HCN_out;
 strout(8).yNHN3 = y_NH3_out;
 strout(8).yCH4 = y_CH4_out;
 strout(8).yH2 = y_H2_out;
 strout(8).yN2 = y_N2_out;

 %x_HCN_in = strin(10).x_HCN; 
x_HCN_in = 0;
 %x_H2O_in = strin(10).x_H2O;
x_H2O_in = 1;
x_in = [x_HCN_in, x_H2O_in];     

x_HCN_out = (G_in_i(2) - y_HCN_out*G_out)/(L_in + G_in_i(2) - y_HCN_out*G_out); 
x_H2O_out = 1-x_HCN_out;
x_out = [x_HCN_out, x_H2O_out];

 strout(9).xHCN = x_HCN_out;
 strout(9).xH2O = x_H2O_out;

G_in_i = G_in.*y_in;
G_out_i= G_out.*y_out;
L_in_i = L_in.*x_in;
L_out_i= L_out.*x_out;

%T_gas_in = 273 + 60;
T_gas_in = strin(7).T 
T_liquid_in = 273 + 20;
t_gas_in = T_gas_in/1000;            % Modified temperature used for the Shomate equations
t_liquid_in = T_liquid_in/1000;

% Standard formation enthalpy (from literature)
 h0f_HCN = cmpin(6).deltaHf ;
 h0f_H2O = cmpin(1).deltaHf ;
 h0f_H2 =  cmpin(5).deltaHf ;
 h0f_N2 =  cmpin(2).deltaHf ;
 h0f_NH3 = cmpin(4).deltaHf ;
 h0f_CH4 = cmpin(3).deltaHf ;

% Heat capacities at 298.15 K (from literature)
 cp_0_HCN = cmpin(6).cp ;      
 cp_0_H2O = cmpin(1).cp ;
 cp_0_H2 =  cmpin(5).cp ;
 cp_0_N2 = cmpin(2).cp;
 cp_0_NH3 = cmpin(4).cp;
 cp_0_CH4 = cmpin(3).cp ;
 
 % Modified formation enthalpies using hf = h0f + cp (T1-T0) with the assumption that the heat capacity is constant

 hf_HCN_gas_in =  h0f_HCN + (cp_0_HCN * (T_gas_in - 293.15));  
 hf_H2_gas_in =  h0f_H2 + (cp_0_H2 * (T_gas_in -  293.15)); 
 hf_N2_gas_in =  h0f_N2 + (cp_0_N2 *(T_gas_in -  293.15)); 
 hf_NH3_gas_in =  h0f_NH3 + (cp_0_NH3*  (T_gas_in -  293.15)); 
 hf_CH4_gas_in = h0f_CH4 + (cp_0_CH4 *(T_gas_in -  293.15)); 

 hf_H2O_liquid_in = h0f_H2O + (cp_0_H2O *(T_liquid_in -  293.15));

 hf_HCN_gas_out = @(T_out) h0f_HCN + (cp_0_HCN * (T_out -  293.15));  
 hf_H2_gas_out = @(T_out) h0f_H2 + (cp_0_H2 * (T_out -  293.15)); 
 hf_N2_gas_out = @(T_out) h0f_N2 + (cp_0_N2 * (T_out - 293.15)); 
 hf_NH3_gas_out = @(T_out) h0f_NH3 + (cp_0_NH3 * (T_out-  293.15)); 
 hf_CH4_gas_out = @(T_out) h0f_CH4 + (cp_0_CH4 * (T_out-  293.15)); 

 hf_HCN_liquid_out = @(T_out) h0f_HCN + (cp_0_HCN * (T_out-  293.15));  
 hf_H20_liquid_out = @(T_out) h0f_H2O + (cp_0_H2O * (T_out -  293.15)); 

% Energy balance with the assumption that the heat capacities are constant

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
Temp_out = fsolve (Ebalance, 300)
Temp_out_celsius = Temp_out - 273
strout(9).T = Temp_out;
strout(8).T = Temp_out;
%%
 % Case 2 : Heat capacitites are not constant anymore
 
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

hf_HCN_in_non = (A(1)*t_gas_in) + (B(1)*((t_gas_in)^2)/2) + (C(1)*((t_gas_in)^3)/3) + (D(1)*((t_gas_in)^4)/4) - (E(1)/t_gas_in) +F(1) - H(1);                      
hf_NH3_in_non = (A(2)*t_gas_in) + (B(2)*((t_gas_in)^2)/2) + (C(2)*((t_gas_in)^3)/3) + (D(2)*((t_gas_in)^4)/4) - (E(2)/t_gas_in) +F(2) - H(2);
hf_H2_in_non = (A(3)*t_gas_in) + (B(3)*((t_gas_in)^2)/2) + (C(3)*((t_gas_in)^3)/3) + (D(3)*((t_gas_in)^4)/4) - (E(3)/t_gas_in) +F(3) - H(3); 
hf_CH4_in_non = (A(4)*t_gas_in) + (B(4)*((t_gas_in)^2)/2) + (C(4)*((t_gas_in)^3)/3) + (D(4)*((t_gas_in)^4)/4) - (E(4)/t_gas_in) +F(4) - H(4); 
hf_N2_in_non = (A(5)*t_gas_in) + (B(5)*((t_gas_in)^2)/2) + (C(5)*((t_gas_in)^3)/3) + (D(5)*((t_gas_in)^4)/4) - (E(5)/t_gas_in) +F(5) - H(5); 

% Enthalpies for the outcoming gas stream

hf_HCN_out_non = @(t_out) (A(1)*t_out) + (B(1)*((t_out)^2)/2) + (C(1)*((t_out)^3)/3) + (D(1)*((t_out)^4)/4) - (E(1)/t_out) +F(1) - H(1);  
hf_NH3_out_non = @(t_out) (A(2)*t_out) + (B(2)*((t_out)^2)/2) + (C(2)*((t_out)^3)/3) + (D(2)*((t_out)^4)/4) - (E(2)/t_out) +F(2) - H(2); 
hf_H2_out_non = @(t_out) (A(3)*t_out) + (B(3)*((t_out)^2)/2) + (C(3)*((t_out)^3)/3) + (D(3)*((t_out)^4)/4) - (E(3)/t_out) +F(3) - H(3);  
hf_CH4_out_non = @(t_out) (A(4)*t_out) + (B(4)*((t_out)^2)/2) + (C(4)*((t_out)^3)/3) + (D(4)*((t_out)^4)/4) - (E(4)/t_out) +F(4) - H(4);  
hf_N2_out_non = @(t_out) (A(5)*t_out) + (B(5)*((t_out)^2)/2) + (C(5)*((t_out)^3)/3) + (D(5)*((t_out)^4)/4) - (E(5)/t_out) +F(5) - H(5);  

% Enthalpies for the outcoming liquid stream

hf_H2O_out_non = @(t_out) (A(7)*t_out) + (B(7)*((t_out)^2)/2) + (C(7)*((t_out)^3)/3) + (D(7)*((t_out)^4)/4) - (E(7)/t_out) +F(7) - H(7);  
hf_HCN_out_non = @(t_out) (A(1)*t_out) + (B(1)*((t_out)^2)/2) + (C(1)*((t_out)^3)/3) + (D(1)*((t_out)^4)/4) - (E(1)/t_out) +F(1) - H(1); 

% Energy balance with the assumption that the heat capacities are constant

 E_gas_in_non = G_in * ((y_HCN_in * hf_HCN_in_non) + (y_H2_in * hf_H2_in_non) + (y_N2_in * hf_N2_in_non) + (y_NH3_in * hf_NH3_in_non) + (y_CH4_in * hf_CH4_in_non));  
 E_liquid_in_non = L_in * (x_H2O_in * hf_H2O_liquid_in);
 E_gas_out_non = @(t_out) G_out * ((y_HCN_out * hf_HCN_out_non(t_out)) + (y_H2_out * hf_H2_out_non(t_out)) + (y_N2_out * hf_N2_out_non(t_out)) + (y_NH3_out * hf_NH3_out_non(t_out)) + (y_CH4_out * hf_CH4_out_non(t_out)));
 E_liquid_out_non = @(t_out) L_out * ((x_H2O_out * hf_H2O_out_non(t_out)) + (x_HCN_out * hf_HCN_out_non(t_out))); 
 
 E_in_non = E_gas_in_non + E_liquid_in_non;
 E_out_non = @(t_out) E_gas_out_non(t_out) + E_liquid_out_non(t_out);
 
 Ebalance_non = @(t_out) (E_out_non(t_out) - E_in_non);
 temp_out_non = fsolve (@(t_out) Ebalance_non(t_out), 298);
 Temp_out_non = temp_out_non * 1000;
 flow_ratio = L_in/G_in;
 flow_ratio_out = L_out/G_out;
 %%
 % Calculation of the height of the HCN absorber column
 % Given values
 
 p = 101325;    % Pressure inside the column in Pa
 M_H2O = 18;
 M_HCN = 27;
 M_H2 = 2;
 M_N2 = 28;
 M_NH3 = 17;
 M_CH4 = 16;
 
 
 % Assumption that the gas phase consists only of hydrogen and HCN
 henry_HCN = HenrysConstant(Temp_out,cmp,6) ;    % Henry coefficient of HCN in [Pa]
 velocity = 2;  % Gas velocity in [m/s]
 dia = diameter_calc2 (velocity,G_in, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);    % Diameter in m
 area = pi*((dia/2)^2);           % Cross sectional area 
 P_bar = 1;                       % Pressure in bar
 mu_sum = mu_calc2(Temp_out);     % Dynamic viscosity [Pa s]
 G_m = G_in*(y_NH3_in*M_NH3 + y_H2_in*M_H2 + y_HCN_in*M_HCN + y_N2_in*M_N2 + y_CH4_in*M_CH4)/1000;
 rho_sum = rho_calc(Temp_out, y_HCN_in, y_H2_in, y_NH3_in, y_CH4_in, y_N2_in)/1000;   % Calculation of the densitites in dependance of the outlet temperature in [g/m^3]
 nu = mu_sum/rho_sum;                                        % Kinematic viscosity [m^2/s]
 dp = 0.025;                                                 % Diameter of the Füllkörper in m
 R = 8.314;                                                  % Ideal Gas constant
 phi = 2.6;                                                  % Empirical parameter for water for the calculation of the diffusioncoefficient
 mol_vol_HCN = (M_HCN/cmpin(6).rho)* 1000;                     % Molar volume of HCN in cm^3/mol
 D_HCN = ((7.4 * 10^(-8)) * Temp_out * ((M_H2O * phi)^2))/(mu_sum * (mol_vol_HCN^0.6)) * 10^(-4);      % Diffusion coefficient in cm^2/s
 Sc = mu_sum/(rho_sum*D_HCN);     % Schmidtzahl 
 v_frei = 0.74;               % Freier Volumenanteil
 tot_surf = 190;              % Total surface in m^2/m^3
 ap = 492;                    % Packungsfaktor trocken in m^2/m^3 
 aw = 587;                    % Packungsfaktor nass aus Perry's
 L_m = L_in * M_H2O/1000;          % Mass flux in kg/s
 g = 9.81;                    % Erdbeschleunigung in m/s
mu_H2O = 547* 10^(-6);        % Dynamic viscosity of water at 50 °C
k_L = ((cmpin(1).rho/(mu_H2O*g))^(-1/3)) * 0.0051 * ((L_m/(aw*mu_H2O))^(2/3)) * (Sc^(-1/2)) * ((ap * dp)^0.4); % Mass transport coefficient in the liquid phase
K_L = k_L *rho_sum/(M_H2O/1000);     % Overall mass transfer coefficient
H_L = L_in/(area * K_L * tot_surf);     % H_L Value
m = henry_HCN/p;   % Gleichgewichtskonstante

%%
% HG Value
Q = G_in*R*Temp_out/P_bar;
sum_V_H2 = 7.07;
sum_V_HCN = 24.17  ;
D = 10^(-3)*Temp_out^1.75*(1/M_H2 + 1/M_HCN)^0.5/(((sum_V_H2^(1/3))+ (sum_V_HCN^(1/3)))^2)*10^(-4);
Sc = mu_sum/(rho_sum*D);
KG = 5.23*ap*D/(R*Temp_out) * (G_m/(ap*nu))^0.7 * Sc^(1/3) * (ap*dp)^(-2)*p;
HG = G_in/(area * KG * tot_surf);
HTU = HG + ((m/L_in/G_in)*H_L);
A = L_in/(G_in * m);
alpha = (y_HCN_in-y_HCN_out)/(y_HCN_in -m*x_HCN_in);
NTU=A./(A-1).*log((1-alpha./A)/(1-alpha));
h = HTU * NTU;
ratio = h/dia;
V_column = (dia/2)^2 * h * pi;
 %%
 % Parameters for the Van der Waals equation with a: [bar L^2/mol^2] and b:[L/mol]
 % HCN, H2, N2, CH4, NH3, H2O
 a = [11.29, 0.2453, 1.370, 2.300, 4.225, 5.537];
 b = [0.08806, 0.02651, 0.0387, 0.04301, 0.03713, 0.03049];
%%
% CAPEX 
CAPEX_hcnabs = capex_calc(h, dia);
CAPEX_mil = CAPEX_hcnabs/1000000;

% OPEX
[opex_H2O, opex_H2SO4, opex_wastewater, opex_tot] = opex_calc(G_in, L_in, x_H2O_in, 0, 0, 0, 0);

%L min bestimmen
min_ratio = ((y_HCN_in - y_HCN_out)/((y_HCN_in/m)-x_HCN_in));
L_min = G_in * min_ratio;
L_true = L_min * 1.5;

fprintf('Number of theoretical units: NTU = %g\n', NTU);
fprintf('Height of theoretical units: HTU = %g\n', HTU);
fprintf('Flow rate ratio: L/G = %g\n', flow_ratio);
fprintf('Outlet temperature [°C]: T = %g\n', Temp_out_celsius);
fprintf('Column Height [m]: H = %g\n', h);
fprintf('Column Diameter [m]: D = %g\n', dia );
fprintf('HCN Column CAPEX [Mio. US$]: Capex = %g\n', CAPEX_mil);
fprintf('HCN Column OPEX [Mio. US$]: Opex = %g\n', opex_tot);

untout(3).h = h;
untout(3).capex = CAPEX_mil;
untout(3).opex = opex_tot;
untout(3).ntu = NTU;
untout(3).htu = HTU;
untout(3).V = V_column;


 end
