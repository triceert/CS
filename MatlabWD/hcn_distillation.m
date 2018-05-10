function [cmpout, untout, strout] = hcn_distillation(cmpin, untin, strin)
% Calculation of HCN distillation using the Fenske-Underwood-Gilliland
% method
% since we are at bubble point (Aufgabenstellung: "f�hren Sie den Feed im Siedezustand zu"), there's only L coming
% in (no V) --> need a HX in stream 9 to get to saturated liquid 

%%%%%%%%%%%%%%%%
% cd /Users/Clemens/CS/MatlabWD
%%%%%%%%%%%%%%%%

%% Setting/extracting parameters
R=untin(5).idgc;                                                       % ideal gas constant [J.mol-1.K-1]
feed_L = strin(9).L;                                                   % molar liquid flowrate in feed                                             
feed_V = strin(9).G;                                                   % molar gas flowrate in feed [mol/hr] (= 0 since only liquid)
q = 1;                                                                 % feed quality 
                                                                       % (fraction of feed that is liquid, q=1 since @ bubble point)
pressure = 1e5;                                                        % operating the column at atmospheric pressure
MW_H2O = cmpin(1).MW;                                                  % molecular weight H2O [kg/mol]
MW_HCN = cmpin(6).MW;                                                  % molecular weight HCN [kg/mol]

if untin(1).ideal_real == 0
    thermo_model = 'ideal';                                            % thermodynamic model, can use 'nrtl' or 'ideal'
else
    thermo_model = 'nrtl';                                             % thermodynamic model, can use 'nrtl' or 'ideal'
end

if strcmp(thermo_model, 'ideal')
    feed_L = strin(9).L;                                                   % molar liquid flowrate in feed                                             
    feed_V = strin(9).G;                                                   % molar gas flowrate in feed [mol/hr] (= 0 since only liquid)
else
    feed_L = strin(9).Lreal;                                               % molar liquid flowrate in feed                                             
    feed_V = strin(9).G;     % doesn't need to be changed since = 0    % molar gas flowrate in feed [mol/hr] (= 0 since only liquid)
end 


T_boiling_H2O = cmpin(1).bp; 
T_boiling_HCN = cmpin(6).bp; 

P_sat_H2O = @(temperature) antoine_equation_new(cmpin, temperature, 1); 
P_sat_HCN = @(temperature) antoine_equation_new(cmpin, temperature, 6);
       



z.HCN = (strin(9).xHCN*feed_L + strin(9).yHCN *feed_V)/(feed_L+feed_V);  
z.H2O = (strin(9).xH2O*feed_L + strin(9).yH2O *feed_V)/(feed_L+feed_V); 
%z.HCN = 0.025; 
%z.H2O = 1-z.HCN; 
%z.AS = (strin(9).xAS*feed_L + strin(9).yAS *feed_V)/(feed_L+feed_V);  
%z.H2 = (strin(9).xH2*feed_L + strin(9).yH2 *feed_V)/(feed_L+feed_V); 
F = feed_L + feed_V; % since the whole outlet stream of the HCN absorption unit will be condensed before entering the column
%F = 531.4; 

%% Specifications for outlet streams
x_LK_D = 0.995;         % from specifications on tasksheet
x_HK_D = 1-x_LK_D; 
x_LK_B = 10e-6; 
x_HK_B = 1-x_LK_B;


%% Calculation of feed conditions & HX before distillation column 
T_feed_before_HX = strin(9).T; 
T_feed = bubblepoint_new([z.H2O, z.HCN], pressure, cmpin, untin, thermo_model); 
[HT_L_before_HX] = enthalpy_temperature_liquid(T_feed_before_HX,cmpin,untin); 
[HT_L_after_HX] = enthalpy_temperature_liquid(T_feed,cmpin,untin);
Q_HX_before_distillation_column = z.H2O*HT_L_after_HX(1)+z.HCN*HT_L_after_HX(6)-z.H2O*HT_L_before_HX(1)-z.HCN*HT_L_before_HX(6); 
T_steam_at_6_bar = 158.8+273.15;                        % [K], 6 bar: see task sheet,                                                        
                                                % from https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
specific_enthalpy_of_evaporation_of_steam_at_6_bar = 2257e3; % [J/kg], 
                                           % from https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
area_HX_before_distillation_column = Q_HX_before_distillation_column/(700*(T_steam_at_6_bar-T_feed_before_HX));
steam_flow_HX_before_distillation_column = Q_HX_before_distillation_column/specific_enthalpy_of_evaporation_of_steam_at_6_bar; % [kg/s]
OPEX_HX_before_distillation_column = steam_flow_HX_before_distillation_column/1000*8000*3600*20; % steam cost: 20 USD/tonne
CAPEX_HX_before_distillation_column = 25000*(area_HX_before_distillation_column^0.65);

% want to operate the column at atmospheric pressure 
% pressure = @(temperature) z.CH4*P_sat_CH4(temperature) + z.NH3*P_sat_NH3(temperature) + z.Egas*P_sat_Egas(temperature) + 
% z.H2SO4*P_sat_H2SO4(temperature) + z.HCN*P_sat_HCN(temperature) + z.AS*P_sat_AS(temperature) + z.H2*P_sat_H2(temperature);
% formula from solution Separation Technologies, Series 6
%pressure = @(temperature) z.HCN*P_sat_HCN(temperature) + z.H2O*P_sat_H2O(temperature);
pressure = strin(9).p; % Pa

% NEED TO START AGAIN HERE FOR PRESSURE DROP CALCULATIONS
if strcmp(thermo_model, 'ideal') 
    alpha = @(temperature, x_LK) P_sat_HCN(temperature)/P_sat_H2O(temperature); % relative volatility
end
if strcmp(thermo_model, 'nrtl') || strcmp(thermo_model, 'vanlaar')
    alpha = @(temperature, x_LK) relative_volatility(x_LK, temperature, cmpin, thermo_model); % relative volatility
end
% --> one sees that alpha doesn't change dramatically
% calculation of geometric mean of alpha: 


alpha_mean = (alpha(T_boiling_H2O, x_LK_B) * alpha(T_boiling_HCN, x_LK_D) * alpha(T_feed, z.HCN))^(1/3);
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
%N_S_real = ceil((Y+N_S_min)/(1-Y)); % rounding up to nearest integer
N_S_real = (Y+N_S_min)/(1-Y);

N_S_real_new = 0; % initializing so that we do at least one iteration

while N_S_real ~= N_S_real_new
    pressure_drop = 700*N_S_real; % 0.7 kPa/tray pressure drop over bubble cap column
                                  % from: "Separation Process Principles", J.
                                  % D. Seader, E. J. Henley, p. 375
    pressure_top = pressure-pressure_drop; 
    new_temperature_at_top = bubblepoint_new([x_HK_D, x_LK_D], pressure_top, cmpin, untin, thermo_model);
    new_alpha_mean = (alpha(T_boiling_H2O, x_LK_B) * alpha(new_temperature_at_top, x_LK_D) * alpha(T_feed, z.HCN))^(1/3);
    N_S_min = log(x_LK_D/x_HK_D*x_HK_B/x_LK_B)/log(new_alpha_mean)-1; 
    %new_theta_sol = fsolve (@(theta) (new_alpha_mean*z.HCN)/(new_alpha_mean-theta)+(1*z.H2O)/(1-theta), theta0, options);
    new_theta_sol = new_alpha_mean*(z.HCN+z.H2O)/(new_alpha_mean*z.HCN+z.H2O);
    RR_min = (new_alpha_mean*x_LK_D)/(new_alpha_mean-new_theta_sol)+(1*x_HK_D)/(1-new_theta_sol)-1; 
    RR_real = RR_min * 1.5; % heuristic, from Stavros' introduction
    X = (RR_real-RR_min)/(RR_real+1); 
    Y = 1-exp(((1+54.4*X)/(11+117.2*X))*((X-1)/sqrt(X)));
    N_S_real = N_S_real_new;                % saving for next iteration   
    N_S_real_new = round((Y+N_S_min)/(1-Y),1);% rounding up to nearest integer
    
    clear new_theta_sol;  
end
% when this while-loop terminates, N_S_real = N_S_real_new --> no need to replace it
efficiency = 0.75;                                      % Stavros said a plate efficiency between 0.7 and 0.8 is reasonable
N_S_real = ceil(N_S_real/efficiency);                   % accounting for plate efficiency <1 and rounding up to nearest integer

 



%%% SIZING THE COLUMN %%% 
height = 1.2*N_S_real*0.6;                              % 0.5 m distance between plates, so use 0.6 m/plate
D = z.HCN*F/0.995;                                      % neglecting the 10 ppm HCN in bottom stream
V_R = (RR_real+1)*D; 
V_S = V_R;                                              % since at q=1
V_flowrate = V_S*R*T_boiling_H2O/pressure;              % using V_S since this refers to stripping section (bottom of column), 
                                                        % which is hottest --> lowest gas density at same pressure
%rho_H2O_B = MW_H2O*pressure/(8.3144*T_boiling_H2O); 
%rho_HCN_B = MW_HCN*pressure/(8.3144*T_boiling_H2O); 
%rho_V = x_LK_B*rho_HCN_B+x_HK_B*rho_H2O_B;             % weighted average of densities
rho_V = 0.598;                                          % in reboiler: almost pure water, 
                                                        % from https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
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
T_cooling_water_in=5;                             % assumed temperature difference of cooling water, [K], goes from 5�C --> 15 �C 
T_cooling_water_out=15;                                 % assumption
delta_T1 = T_boiling_HCN - T_cooling_water_in;          % needed for LMTD calculation
delta_T2 = T_boiling_HCN - T_cooling_water_out;         % needed for LMTD calculation
deltaH_HCN = cmpin(6).Hv;                               % Molar enthalpy of vaporization of HCN
L_0 = D*RR_real; 
V_1 = L_0 + D; 
Q_cond = deltaH_HCN*V_1;                                % Heat duty condenser [J/s], "-" due to convention
LMTD = (delta_T1-delta_T2)/log(delta_T1/delta_T2);      % Using log-mean temperature difference for counter-current HX
area_cond = Q_cond/(700*LMTD);                          % 0.700 kW/(m2*K) from task sheet
enthalpy_feed = enthalpy_temperature_liquid(T_feed,cmpin,untin); 
enthalpy_distillate = enthalpy_temperature_liquid(new_temperature_at_top,cmpin,untin); 
enthalpy_bottom = enthalpy_temperature_liquid(T_boiling_H2O,cmpin,untin); 
h_F = z.H2O*enthalpy_feed(1)+z.HCN*enthalpy_feed(6); 
h_D = x_HK_D*enthalpy_distillate(1) + x_LK_D*enthalpy_distillate(6);
h_B = x_HK_B*enthalpy_bottom(1) + x_LK_B*enthalpy_bottom(6);        % not yet any residual enthalpies
Q_reboiler = D*h_D+B*h_B-F*h_F+Q_cond;               % "-" due to convention
T_steam_at_6_bar = 158.8+273.15;                        % [K], 6 bar: see task sheet,                                                        
                                                % from https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
area_reboiler = Q_reboiler/(700*(T_steam_at_6_bar-T_boiling_H2O));
heat_capacity_cooling_water = heat_capacity((new_temperature_at_top+T_cooling_water_in)/2,cmpin,untin, 1); 
cP_cooling_water = @(temperature) heat_capacity(temperature, cp_coefficients_cooling_water);
cooling_water_mass_flow = Q_cond*MW_H2O/(heat_capacity_cooling_water*(new_temperature_at_top-T_cooling_water_in)); 

% cooling water mass flow [kg/s], evaluating cp_cooling_water at average temperature (cp of H2O doesn't change much in the 
% region in question anyways) 
specific_enthalpy_of_evaporation_of_steam_at_6_bar = 2257e3; % [J/kg], 
                                           % from https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
steam_flow_condenser = Q_reboiler/(specific_enthalpy_of_evaporation_of_steam_at_6_bar); % [kg/s] steam mass flow rate



%%% Cost calculation
OPEX_cooling_water = (cooling_water_mass_flow/1000)*8000*3600*0.15; 
% using 0.15 USD/tonne instead of the given 0.10 USD/tonne to account
                                                               % for the need to pre-cool
OPEX_steam = steam_flow_condenser/1000*8000*3600*20;           % steam cost: 20 USD/tonne
                                                        % [USD/a], /1000 to convert to tonnes, 
                                                        % *8000 since 8000 operating hours/a
CAPEX_column_itself = 80320*(height^0.76)*(d_min_bottom^1.21);
CAPEX_cond = 25000*(area_cond^0.65); 
CAPEX_reboiler = 25000*(area_reboiler^0.65); 
%OPEX_column = OPEX_cooling_water + OPEX_steam; 
%%% WHICH TEMPERATURE DIFFERENCE FOR REBOILER HEATED WITH STEAM? 



%% OUTPUTS
cmpout = cmpin; 
untout = untin; 
strout = strin; 

untout(4).h = height; 
untout(4).rad = d_min_bottom/2; 
untout(4).V = (pi*(untout(4).rad)^2)*untout(4).h; 
untout(4).En = Q_cond + Q_reboiler; 
untout(4).capex = CAPEX_column_itself+CAPEX_cond+CAPEX_reboiler+CAPEX_HX_before_distillation_column; 
untout(4).opex = OPEX_steam + OPEX_cooling_water+OPEX_HX_before_distillation_column; 
strout(11).L = D; 
strout(11).G = 0; % total condenser
strout(11).p = pressure_top; 
strout(11).T = new_temperature_at_top; 
strout(11).xHCN=x_LK_D;
strout(11).xH2O=x_HK_D;
strout(11).xCH4=0;
strout(11).xNH3=0;
strout(11).xEgas=0;
strout(11).xAS=0;
strout(11).xH2=0;
strout(11).xN2=0;
strout(11).xH2SO4=0;  
strout(11).yHCN=0;     % since Stream 11 is all liquid
strout(11).yH2O=0;     % since Stream 11 is all liquid
strout(11).yCH4=0;     % since Stream 11 is all liquid
strout(11).yNH3=0;     % since Stream 11 is all liquid
strout(11).yEgas=0;    % since Stream 11 is all liquid
strout(11).yAS=0;      % since Stream 11 is all liquid
strout(11).yH2=0;      % since Stream 11 is all liquid
strout(11).yN2=0;      % since Stream 11 is all liquid
strout(11).yH2SO4=0;   % since Stream 11 is all liquid
if strcmp(thermo_model, 'ideal')
    strout(10).L = B; 
else
    strout(10).Lreal = B; 
end
strout(10).G = 0;      % outlet stream 10 is all liquid
strout(10).p = pressure; 
strout(10).T = T_boiling_H2O; 
strout(10).xHCN=x_LK_B;
strout(10).xH2O=x_HK_B;
strout(10).xCH4=0;
strout(10).xNH3=0;
strout(10).xEgas=0;
strout(10).xAS=0;
strout(10).xH2=0;
strout(10).xN2=0;
strout(10).xH2SO4=0;  
strout(10).yHCN=0;     % since Stream 11 is all liquid
strout(10).yH2O=0;     % since Stream 11 is all liquid
strout(10).yCH4=0;     % since Stream 11 is all liquid
strout(10).yNH3=0;     % since Stream 11 is all liquid
strout(10).yEgas=0;    % since Stream 11 is all liquid
strout(10).yAS=0;      % since Stream 11 is all liquid
strout(10).yH2=0;      % since Stream 11 is all liquid
strout(10).yN2=0;      % since Stream 11 is all liquid
strout(10).yH2SO4=0;   % since Stream 11 is all liquid






if 0
%%% Plots
[~, ~, raw] = xlsread('Txy-HCN-H2O.xlsx','Txy HCN-H2O');
raw = raw(2:end,:);
data = reshape([raw{:}],size(raw));
plot_x_HCN = data(:,1);
plot_y_HCN = data(:,2);
plot_T = data(:,3);
clearvars data raw;
figure(3); 
scatter(plot_x_HCN, plot_T-273); 
hold on; 
scatter(plot_y_HCN, plot_T-273); 

x1=0:0.01:1; 
for i=1:101
    y1(i) = bubblepoint_new([x1(i), 1-x1(i)], 1e5, cmpin, untin, thermo_model)-273; 
end
plot(1-x1, y1);
legend('Bubble point experimental', 'Dew point experimental', 'Bubble point model');
end
%N_S_min
%N_S_real
%RR_min
%RR_real

%test = 1;



end

