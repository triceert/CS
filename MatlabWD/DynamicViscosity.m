function [pmu] = DynamicViscosity(cmp)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: cmp = compound struct     
%
% OUTPUT: plambda 9x3 array with polynom of the compound on each line 
% Indexes: 1: Water, 2:Nitrogen,3:Methane, 4:Ammonia,5: Hydrogen, 6:Hydrogen cyanide,
% 7:Sulfuric acid, 8:Ammonium sulfate, 9:Carbon dioxide

%Gases
Tgas100 = 100:100:600; % [K]
Tgas200 = 200:100:600; % [K]
Tgas300 = 300:100:600; % [K]
Tgas400 = 400:100:600; % [K]

mugas1 = [cmp(1).mugas400,cmp(1).mugas500,cmp(1).mugas600];
mugas2 = [cmp(2).mugas100,cmp(2).mugas200,cmp(2).mugas300,cmp(2).mugas400,cmp(2).mugas500,cmp(2).mugas600];
mugas3 = [cmp(3).mugas100,cmp(3).mugas200,cmp(3).mugas300,cmp(3).mugas400,cmp(3).mugas500,cmp(3).mugas600];
mugas4 = [cmp(4).mugas300,cmp(4).mugas400,cmp(4).mugas500,cmp(4).mugas600];
mugas5 = [cmp(5).mugas100,cmp(5).mugas200,cmp(5).mugas300,cmp(5).mugas400,cmp(5).mugas500,cmp(5).mugas600];
mugas9 = [cmp(9).mugas200,cmp(9).mugas300,cmp(9).mugas400,cmp(9).mugas500,cmp(9).mugas600];

%Liquids
%Tliq1 = 273.16:25:373.16; At those temperatures, water will not be liquid
Tliq6 = [273.16,298.16];
Tliq7 = 298.16:25:373.16;

%muliq1 = [cmp(1).muliq273,cmp(1).muliq298,cmp(1).muliq323,cmp(1).muliq348,cmp(1).muliq373];
muliq6 = [cmp(6).muliq273,cmp(6).muliq298];
muliq7 = [cmp(7).muliq298,cmp(7).muliq323,cmp(7).muliq348,cmp(7).muliq373];

%Polynomial fitting
pmu = zeros(9,3);

%pmu(1,:) = polyfit([Tliq1 Tgas400],[muliq1 mugas1],2);
pmu(1,:) = polyfit(Tgas400,mugas1,2);
pmu(2,:) = polyfit(Tgas100,mugas2,2);
pmu(3,:) = polyfit(Tgas100,mugas3,2);
pmu(4,:) = polyfit(Tgas300,mugas4,2);
pmu(5,:) = polyfit(Tgas100,mugas5,2);
pmu(6,:) = polyfit(Tliq6,muliq6,2);
pmu(7,:) = polyfit(Tliq7,muliq7,2);
pmu(8,:) = [NaN NaN NaN];
pmu(9,:) = polyfit(Tgas200,mugas9,2);
end