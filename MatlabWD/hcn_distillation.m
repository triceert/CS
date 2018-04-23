function [cmpout, untout, strout] = hcn_distillation(cmpin, untin, strin)
% Calculation of HCN distillation using the Fenske-Underwood-Gilliland
% method
% since we are at bubble point (Aufgabenstellung: "führen Sie den Feed im Siedezustand zu"), there's only L coming
% in (no V) --> might need a HX in stream 9 to get to saturated liquid 

%%%%%%%%%%%%%%%%
% cd /Users/Clemens/CS/MatlabWD
%%%%%%%%%%%%%%%%

feed_L = strin(9).L 
feed_V = strin(9).G; % molar flowrate in feed [mol/hr]

q = 1; % feed quality (fraction of feed that is liquid, q=1 since @ bubble point)
z.H2O = (strin(9).xH2O*feed_L + strin(9).yH2O *feed_V)/(feed_L+feed_V); 
z.CH4 = (strin(9).xCH4*feed_L + strin(9).yCH4 *feed_V)/(feed_L+feed_V); 
z.NH3 = (strin(9).xNH3*feed_L + strin(9).yNH3 *feed_V)/(feed_L+feed_V);  
z.Egas = (strin(9).xEgas*feed_L + strin(9).yEgas *feed_V)/(feed_L+feed_V);  
z.H2SO4 = (strin(9).xH2SO4*feed_L + strin(9).yH2SO4 *feed_V)/(feed_L+feed_V);  
z.HCN = (strin(9).xHCN*feed_L + strin(9).yHCN *feed_V)/(feed_L+feed_V); 
z.AS = (strin(9).xAS*feed_L + strin(9).yAS *feed_V)/(feed_L+feed_V);  
z.H2 = (strin(9).xH2*feed_L + strin(9).yH2 *feed_V)/(feed_L+feed_V); 
F = feed_L + feed_V; % since the whole outlet stream of the HCN absorption unit will be condensed before entering the column


A_CH4 = 6.7021; 	% Antoine parameters from http://ddbonline.ddbst.com/AntoineCalculation/AntoineCalculationCGI.exe?component=Methane 
B_CH4 = 394.48; 	 
C_CH4 = 264.609;

A_NH3 = 4.86886;   % https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
B_NH3 = 1113.928;
C_NH3 = -10.409; 

A_H2O = 7.96681; B_H2O = 1668.21; C_H2O = 228.0; % (values for T > 60 °C), from Excel file  
A_HCN = 7.5282; B_HCN = 1329.5; C_HCN = 260.4;      % from Excel file

P_sat_H2O = @(temperature) antoine_equation(A_H2O, B_H2O, C_H2O, temperature);
P_sat_CH4 = @(temperature) antoine_equation(A_CH4, B_CH4, C_CH4, temperature); 
P_sat_NH3 = @(temperature) antoine_equation(A_NH3, B_NH3, C_NH3, temperature); 
P_sat_Egas = @(temperature) 10^10; 
P_sat_H2SO4 = @(temperature) 10^10;
P_sat_HCN = @(temperature) antoine_equation(A_HCN, B_HCN, C_HCN, temperature);
P_sat_AS = @(temperature) 10^10;
P_sat_H2 = @(temperature) 10^10;

T_boiling_H2O = 373.15; 
T_boiling_HCN = 299; % from Wikipedia, [K]
T_boiling_mixture_assumption = 373.15*0.8+299*0.2; % just a shitty weighted average, no thermodynamics here 
%z.CH4 = 0.4; z.NH3 = 0.1; z.Egas = 0.1; z.H2SO4 = 0.1; z.HCN = 0.1; z.AS = 0.1; z.H2 = 0.1; 
z.H2O = 0.8; z.HCN = 0.2; 
% want to operate the column at atmospheric pressure 
% pressure = @(temperature) z.CH4*P_sat_CH4(temperature) + z.NH3*P_sat_NH3(temperature) + z.Egas*P_sat_Egas(temperature) + z.H2SO4*P_sat_H2SO4(temperature) + z.HCN*P_sat_HCN(temperature) + z.AS*P_sat_AS(temperature) + z.H2*P_sat_H2(temperature);
% formula from solution Separation Technologies, Series 6
pressure = @(temperature) z.HCN*P_sat_HCN(temperature) + z.H2O*P_sat_H2O(temperature);
alpha = @(temperature) P_sat_HCN(temperature)/P_sat_H2O(temperature); % relative volatility
% --> one sees that alpha doesn't change dramatically
% calculation of geometric mean of alpha: 
alpha_mean = (alpha(T_boiling_H2O) * alpha(T_boiling_HCN) * alpha(T_boiling_mixture_assumption))^(1/3);

%%% Fenske equation %%%
x_LK_D = 0.995; 
x_HK_D = 1-x_LK_D; 
x_LK_B = 10e-6; 
x_HK_B = 1-x_LK_B; 

N_S_min = log(x_LK_D/x_HK_D*x_HK_B/x_LK_B)/log(alpha_mean)-1; 

%%% Underwood formula %%% 
theta0 = 1.1;
theta_sol = fsolve (@(theta) (alpha_mean*z.HCN)/(alpha_mean-theta)+(1*z.H2O)/(1-theta), theta0);
% theta_sol = 1.6183, so between 1 and alpha_mean --> :) 
RR_min = (alpha_mean*x_LK_D)/(alpha_mean-theta_sol)+(1*x_HK_D)/(1-theta_sol)-1; 
RR_real = RR_min * 1.5; % heuristic, from Stavros' introduction


%%% Gilliland correlation %%%
X = (RR_real-RR_min)/(RR_real+1); 
Y = 1-exp(((1+54.4*X)/(11+117.2*X))*((X-1)/sqrt(X))); 
N_S_real = (Y+N_S_min)/(1-Y); 



end

