function [cmpout, untout, strout] = hcn_distillation(cmpin, untin, strin, thermo_model)
% Calculation of HCN distillation using the Fenske-Underwood-Gilliland
% method
% since we are at bubble point (Aufgabenstellung: "f�hren Sie den Feed im Siedezustand zu"), there's only L coming
% in (no V) --> might need a HX in stream 9 to get to saturated liquid 

%%%%%%%%%%%%%%%%
% cd /Users/Clemens/CS/MatlabWD
%%%%%%%%%%%%%%%%
R=untin(5).idgc;                                                           % ideal gas constant [J.mol-1.K-1]
feed_L = strin(9).L;                                                    
feed_V = strin(9).G;                                                       % molar flowrate in feed [mol/hr]
q = 1;                                                                     % feed quality 
pressure = 1e5; 
                                                                           
                                                                           % (fraction of feed that is liquid, q=1 since @ bubble point)
z.H2O = (strin(9).xH2O*feed_L + strin(9).yH2O *feed_V)/(feed_L+feed_V); 
z.HCN = (strin(9).xHCN*feed_L + strin(9).yHCN *feed_V)/(feed_L+feed_V); 
%z.HCN = 0.035; 
%z.H2O = 1-z.HCN; 
%z.AS = (strin(9).xAS*feed_L + strin(9).yAS *feed_V)/(feed_L+feed_V);  
%z.H2 = (strin(9).xH2*feed_L + strin(9).yH2 *feed_V)/(feed_L+feed_V); 
F = feed_L + feed_V; % since the whole outlet stream of the HCN absorption unit will be condensed before entering the column
F = 100; 
MW_H2O = cmpin(1).MW; 
MW_HCN = cmpin(6).MW;
thermo_model = 'nrtl';  

A_H2O = cmpin(1).anta100; B_H2O = cmpin(1).antb100; C_H2O = cmpin(1).antc100; 
A_HCN = cmpin(6).anta100; B_HCN = cmpin(6).antb100; C_HCN = cmpin(6).antc100; 

P_sat_H2O = @(temperature) antoine_equation(A_H2O, B_H2O, C_H2O, temperature);
P_sat_HCN = @(temperature) antoine_equation(A_HCN, B_HCN, C_HCN, temperature);

x_LK_D = 0.995;         % from specifications on tasksheet
x_HK_D = 1-x_LK_D; 
x_LK_B = 10e-6; 
x_HK_B = 1-x_LK_B;


T_boiling_H2O = cmpin(1).bp; 
T_boiling_HCN = cmpin(6).bp;  
%T_boiling_mixture_assumption = 373.15*0.8+299*0.2; % just a shitty weighted average, no thermodynamics here 
%--> use feed temperature once available 
T_boiling_mixture_assumption = 365.5; %%% need to do bubble point calculation
T_feed = bubblepoint([z.H2O, z.HCN], pressure, cmpin, untin); 



% want to operate the column at atmospheric pressure 
% pressure = @(temperature) z.CH4*P_sat_CH4(temperature) + z.NH3*P_sat_NH3(temperature) + z.Egas*P_sat_Egas(temperature) + 
% z.H2SO4*P_sat_H2SO4(temperature) + z.HCN*P_sat_HCN(temperature) + z.AS*P_sat_AS(temperature) + z.H2*P_sat_H2(temperature);
% formula from solution Separation Technologies, Series 6
%pressure = @(temperature) z.HCN*P_sat_HCN(temperature) + z.H2O*P_sat_H2O(temperature);
pressure = strin(9).p; % Pa


if strcmp(thermo_model, 'ideal') 
    alpha = @(temperature, x_LK) P_sat_HCN(temperature)/P_sat_H2O(temperature); % relative volatility
end
if strcmp(thermo_model, 'nrtl') 
    alpha = @(temperature, x_LK) relative_volatility(x_LK, temperature); % relative volatility
end
% --> one sees that alpha doesn't change dramatically
% calculation of geometric mean of alpha: 


alpha_mean = (alpha(T_boiling_H2O, x_LK_B) * alpha(T_boiling_HCN, x_LK_D) * alpha(T_boiling_mixture_assumption, z.HCN))^(1/3);
% always need to pass alpha the mole fraction of the light key (--> HCN) at
% the location corresponding to the used temperature

%%% FENSKE EQUATION %%% 
N_S_min = log(x_LK_D/x_HK_D*x_HK_B/x_LK_B)/log(alpha_mean)-1; 

%%% UNDERWOOD FORMULA %%%
options = optimset('Display','off');
theta0 = 1.1;                      
theta_sol = fsolve (@(theta) (alpha_mean*z.HCN)/(alpha_mean-theta)+(1*z.H2O)/(1-theta), theta0, options);
% theta_sol between 1 and alpha_mean --> :) 
RR_min = (alpha_mean*x_LK_D)/(alpha_mean-theta_sol)+(1*x_HK_D)/(1-theta_sol)-1; 
RR_real = RR_min * 1.5; % heuristic, from Stavros' introduction


%%% GILLILAND CORRELATION %%%
X = (RR_real-RR_min)/(RR_real+1); 
Y = 1-exp(((1+54.4*X)/(11+117.2*X))*((X-1)/sqrt(X))); 
N_S_real = (Y+N_S_min)/(1-Y); 


%%% SIZING THE COLUMN %%% 
height = 1.2*N_S_real*0.6;                              % 0.5 m distance between plates, so use 0.6 m/plate
D = z.HCN*F/0.995;                                      % neglecting the 10 ppm HCN in bottom stream
V_R = (RR_real+1)*D; 
V_S = V_R;                                              % since at q=1
V_flowrate = V_S*R*T_boiling_H2O/pressure;              % using V_S since this refers to stripping section (bottom of column), 
                                                        % which is hottest --> lowest gas density at same pressure
rho_H2O_B = MW_H2O*pressure/(8.3144*T_boiling_H2O); 
rho_HCN_B = MW_HCN*pressure/(8.3144*T_boiling_H2O); 
rho_V = x_LK_B*rho_HCN_B+x_HK_B*rho_H2O_B;              % weighted average of densities
rho_V_imperial_units = 0.06242796*rho_V;                % converting density from kg/m3 to lbm/ft3
v_max_imperial_units = 1/sqrt(rho_V_imperial_units);    % empirical formula only works with v_max [ft/s] and rho_V [lbm/ft3]
v_max = 0.3048*v_max_imperial_units;                    % converting from ft/s to m/s
A_o_bottom = V_flowrate/v_max/0.6;                      % cross-sectional area of column [m2]; 
                                                        % 0.6 since only 60% avail. for flow (see task sheet)
d_min_bottom = sqrt(4*A_o_bottom/pi);                   % column diameter [m]

%%% MOLAR FLOWRATES IN STRIPPING AND RECTIFICATION SECTION %%% 
L_R = D*RR_real; 
L_S = F+L_R; 
B = F-D; 

%%% Cost condenser & reboiler 
% neglecting water in distillate: 
T_cooling_water_in=5;                                   % assumed temperature difference of cooling water, [K], goes from 5�C --> 15 �C 
T_cooling_water_out=15;                                 % assumption
delta_T1 = T_boiling_HCN - T_cooling_water_in;          % needed for LMTD calculation
delta_T2 = T_boiling_HCN - T_cooling_water_out;         % needed for LMTD calculation
deltaH_HCN = cmpin(6).Hv;                               % Molar enthalpy of vaporization of HCN
Q_cond = deltaH_HCN*V_R;                                % Heat duty condenser [J/s]
LMTD = (delta_T1-delta_T2)/log(delta_T1/delta_T2);      % Using log-mean temperature difference for counter-current HX
area_cond = Q_cond/(700*LMTD);                          % 0.700 kW/(m2*K) from task sheet
enthalpy_feed = enthalpy_temperature_liquid(T_boiling_mixture_assumption,cmpin,untin); 
enthalpy_distillate = enthalpy_temperature_liquid(T_boiling_HCN,cmpin,untin); 
enthalpy_bottom = enthalpy_temperature_liquid(T_boiling_H2O,cmpin,untin); 
h_F = z.H2O*enthalpy_feed(1)+z.HCN*enthalpy_feed(6); 
h_D = x_HK_D*enthalpy_distillate(1) + x_LK_D*enthalpy_distillate(6);
h_B = x_HK_B*enthalpy_bottom(1) + x_LK_B*enthalpy_bottom(6);
Q_reboiler = D*h_D+B*h_B-F*h_F-Q_cond; 

heat_capacity_cooling_water_all = heat_capacity((T_boiling_HCN+T_cooling_water)/2,cmpin,untin); 
heat_capacity_cooling_water = heat_capacity_cooling_water_all(1); 

cP_cooling_water = @(temperature) heat_capacity(temperature, cp_coefficients_cooling_water);
cooling_water_mass_flow = Q_cond*MW_H2O/(heat_capacity_cooling_water*(T_boiling_HCN-T_cooling_water)); 
% cooling water mass flow [kg/s], evaluating cp_cooling_water at average temperature (cp of H2O doesn't change much in the 
% region in question anyways)



%%% McCabe Thiele
%x_plot = 0:0.01:1; 
%y_diagonal = x_plot; 


%%% Cost calculation
OPEX_cooling_water = (cooling_water_mass_flow/1000)*8000*0.15; % using 0.15 USD/tonne instead of the given 0.10 USD/tonne to account
                                                               % for the need to pre-cool

% [USD/a], /1000 to convert to tonnes, *8000 since 8000 operating hours/a, *0.1 is cost in USD/tonne
CAPEX_column_itself = 80320*(height^0.76)*(d_min_bottom^1.21);
CAPEX_cond = 25000*(area_cond^0.65); 
%OPEX_column = OPEX_cooling_water + OPEX_steam; 
%%% WHICH TEMPERATURE DIFFERENCE FOR REBOILER HEATED WITH STEAM? 



cmpout = cmpin; 
untout = untin; 
strout = strin; 
%untout(4).h = height; 
%untout(4).rad = d_min_bottom/2; 
%untout(4).V = (pi*(untout(4).rad)^2)*untout(4).h; 
%untout(4).En = Q_cond + Q_reboiler; 

test = 1;






end

