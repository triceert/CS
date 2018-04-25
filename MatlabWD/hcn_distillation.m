function [cmpout, untout, strout] = hcn_distillation(cmpin, untin, strin, thermo_model)
% Calculation of HCN distillation using the Fenske-Underwood-Gilliland
% method
% since we are at bubble point (Aufgabenstellung: "führen Sie den Feed im Siedezustand zu"), there's only L coming
% in (no V) --> might need a HX in stream 9 to get to saturated liquid 

%%%%%%%%%%%%%%%%
% cd /Users/Clemens/CS/MatlabWD
%%%%%%%%%%%%%%%%

feed_L = strin(9).L;
feed_V = strin(9).G;                                                    % molar flowrate in feed [mol/hr]
q = 1;                                                                  % feed quality 
% (fraction of feed that is liquid, q=1 since @ bubble point)
z.H2O = (strin(9).xH2O*feed_L + strin(9).yH2O *feed_V)/(feed_L+feed_V); 
z.HCN = (strin(9).xHCN*feed_L + strin(9).yHCN *feed_V)/(feed_L+feed_V); 
%z.AS = (strin(9).xAS*feed_L + strin(9).yAS *feed_V)/(feed_L+feed_V);  
%z.H2 = (strin(9).xH2*feed_L + strin(9).yH2 *feed_V)/(feed_L+feed_V); 
F = feed_L + feed_V; % since the whole outlet stream of the HCN absorption unit will be condensed before entering the column
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
T_boiling_mixture_assumption = 317.17; %%% need to do bubble point calculation


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

%%% Fenske equation %%% 
N_S_min = log(x_LK_D/x_HK_D*x_HK_B/x_LK_B)/log(alpha_mean)-1; 

%%% Underwood formula %%%
options = optimset('Display','off');
theta0 = 1.1; % just some guess
theta_sol = fsolve (@(theta) (alpha_mean*z.HCN)/(alpha_mean-theta)+(1*z.H2O)/(1-theta), theta0, options);
% theta_sol between 1 and alpha_mean --> :) 
RR_min = (alpha_mean*x_LK_D)/(alpha_mean-theta_sol)+(1*x_HK_D)/(1-theta_sol)-1; 
RR_real = RR_min * 1.5; % heuristic, from Stavros' introduction


%%% Gilliland correlation %%%
X = (RR_real-RR_min)/(RR_real+1); 
Y = 1-exp(((1+54.4*X)/(11+117.2*X))*((X-1)/sqrt(X))); 
N_S_real = (Y+N_S_min)/(1-Y); 


%%% Sizing the column %%% 
height = 1.2*N_S_real*0.6;
D = z.HCN*F/0.995; % neglecting the 10 ppm HCN in bottom stream
V_R = (RR_real+1)*D; 
V_S = V_R;                                              % since at q=1
V_flowrate = V_S*8.3144*T_boiling_H2O/pressure;         % using V_S since this refers to stripping section (bottom of column), 
% which is hottest --> lowest gas density at same pressure
rho_H2O_B = MW_H2O*pressure/(8.3144*T_boiling_H2O); 
rho_HCN_B = MW_HCN*pressure/(8.3144*T_boiling_H2O); 
rho_V = x_LK_B*rho_HCN_B+x_HK_B*rho_H2O_B;              % weighted average of densities
rho_V_imperial_units = 0.06242796*rho_V;                % converting density from kg/m3 to lbm/ft3
v_max_imperial_units = 1/sqrt(rho_V_imperial_units);    % empirical formula only works with v_max [ft/s] and rho_V [lbm/ft3]
v_max = 0.3048*v_max_imperial_units;                    % converting from ft/s to m/s
A_o_bottom = V_flowrate/v_max;                          % cross-sectional area of column [m2]
% --> probably need to account for "Anteil der Oberfläche für Volumenfluss:
% 60%" from task sheet
d_min_bottom = sqrt(4*A_o_bottom/pi);                   % column diameter [m]

%%% Get molar flowrates in stripping and rectification section 
L_R = D*RR_real; 
L_S = F+L_R; 
B = F-D; 


%LV_0 = [100 100 100 100];
%LV = fsolve(@(LV) [(LV(1)-LV(2))/F-1; 
%              (LV(4)-LV(3))/F; 
%              (RR_real+1)*D-LV(4);
%              LV(2)/D-RR_real], LV_0, options); 
% LV(1) = L_S, LV(2) = L_R, LV(3) = V_S, LV(4) = V_R

%%% NRTL (from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1))
%delta_g12 = 1298.9610; 
%delta_g21 = 539.9577; 
%alpha12 = 0.3836; 
% need to call nrtl.m with species1 = HCN, species2 = water (see source)
%gamma = nrtl(0.8, T_boiling_mixture_assumption, delta_g12, delta_g21, alpha12); 
%gamma_HCN = gamma(1); 
%gamma_H2O = gamma(2); 




%%% Cost condenser & reboiler 
% neglecting water in distillate: 
deltaH_HCN = cmpin(6).Hv;
Q_cond = deltaH_HCN*D; % Heat duty condenser [J/s]
Q_reboiler = 0; % NEED TO BE FIXED

T_cooling_water=15+273.15; % assumed temperature of available cooling water, [K]
heat_capacity_cooling_water_all = heat_capacity((T_boiling_HCN+T_cooling_water)/2,cmpin,untin); 
heat_capacity_cooling_water = heat_capacity_cooling_water_all(1); 

cP_cooling_water = @(temperature) heat_capacity(temperature, cp_coefficients_cooling_water);
cooling_water_mass_flow = Q_cond*MW_H2O/(heat_capacity_cooling_water*(T_boiling_HCN-T_cooling_water)); 
% cooling water mass flow [kg/s], evaluating cp_cooling_water at average temperature (cp of H2O doesn't change much in the 
% region in question anyways)

test = 1;

%%% McCabe Thiele
x_plot = 0:0.01:1; 
y_diagonal = x_plot; 


%%% Cost calculation
OPEX_cooling_water = (cooling_water_mass_flow/1000)*8000*0.1; 
% [USD/a], /1000 to convert to tonnes, *8000 since 8000 operating hours/a, *0.1 is cost in USD/tonne
CAPEX_column_itself = 80320*(height^0.76)*(d_min_bottom^1.21);

OPEX_column = OPEX_cooling_water + OPEX_steam; 




cmpout = cmpin; 
untout = untin; 
strout = strin; 
%untout(4).h = height; 
%untout(4).rad = d_min_bottom/2; 
%untout(4).V = (pi*(untout(4).rad)^2)*untout(4).h; 
%untout(4).En = Q_cond + Q_reboiler; 








end

