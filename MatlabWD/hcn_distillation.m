function [cmpout, untout, strout] = hcn_distillation(cmpin, untin, strin)
% Calculation of HCN distillation using the Fenske-Underwood-Gilliland
% method
% since we are at bubble point (Aufgabenstellung: "führen Sie den Feed im Siedezustand zu"), there's only L coming
% in (no V) --> if this conflicts with the outlet of the HCN absorption unit (shouldn't), we would need a condenser 
feed_L = strin(9).L 
feed_V = strin(9).G; % molar flowrate in feed [mol/hr]

q = 1; % since @ bubble point
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

P_sat_CH4 = @(temperature) antoine_equation(A_CH4, B_CH4, C_CH4, temperature); 
P_sat_NH3 = @(temperature) antoine_equation(A_NH3, B_NH3, C_NH3, temperature); 
P_sat_Egas = @(temperature) 10^10; 
P_sat_H2SO4 = @(temperature) 10^10;
P_sat_HCN = @(temperature) 10^10;
P_sat_AS = @(temperature) 10^10;
P_sat_H2 = @(temperature) 10^10;



z.CH4 = 0.4; z.NH3 = 0.1; z.Egas = 0.1; z.H2SO4 = 0.1; z.HCN = 0.1; z.AS = 0.1; z.H2 = 0.1; 
% want to operate the column at atmospheric pressure 
pressure = @(temperature) z.CH4*P_sat_CH4(temperature) + z.NH3*P_sat_NH3(temperature) + z.Egas*P_sat_Egas(temperature) + z.H2SO4*P_sat_H2SO4(temperature) + z.HCN*P_sat_HCN(temperature) + z.AS*P_sat_AS(temperature) + z.H2*P_sat_H2(temperature);
% formula from solution Separation Technologies, Series 6

alpha = P_sat_HCN(temperature)/P_sat_H2O(temperature); % relative volatility








end

