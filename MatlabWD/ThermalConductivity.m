function [plambda] = ThermalConductivity(cmp)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: cmp = compound struct     
%
% OUTPUT: plambda 9x2 array with polynom of the compound on each line 
% Indexes: 1: Water, 2:Nitrogen,3:Methane, 4:Ammonia,5: Hydrogen, 6:Hydrogen cyanide,
% 7:Sulfuric acid, 8:Ammonium sulfate, 9:Carbon dioxide

%Gases
Tgas100 = 100:100:600; % [K]
Tgas200 = 200:100:600; % [K]
Tgas300 = 300:100:600; % [K]
Tgas400 = 400:100:600; % [K]

lambdagas1 = [cmp(1).lambdagas400,cmp(1).lambdagas500,cmp(1).lambdagas600];
lambdagas2 = [cmp(2).lambdagas100,cmp(2).lambdagas200,cmp(2).lambdagas300,cmp(2).lambdagas400,cmp(2).lambdagas500,cmp(2).lambdagas600];
lambdagas3 = [cmp(3).lambdagas100,cmp(3).lambdagas200,cmp(3).lambdagas300,cmp(3).lambdagas400,cmp(3).lambdagas500,cmp(3).lambdagas600];
lambdagas4 = [cmp(4).lambdagas300,cmp(4).lambdagas400,cmp(4).lambdagas500,cmp(4).lambdagas600];
lambdagas5 = [cmp(5).lambdagas100,cmp(5).lambdagas200,cmp(5).lambdagas300,cmp(5).lambdagas400,cmp(5).lambdagas500,cmp(5).lambdagas600];
lambdagas9 = [cmp(9).lambdagas200,cmp(9).lambdagas300,cmp(9).lambdagas400,cmp(9).lambdagas500,cmp(9).lambdagas600];

%Other compounds

%Tliq = 273.16:25:373.16; It will not be liquid
%lambdaliq1 = [cmp(1).lambdaliq273,cmp(1).lambdaliq298,cmp(1).lambdaliq323,cmp(1).lambdaliq348,cmp(1).lambdaliq373];

TLange = 253.16:20:333.16;
lambdaliq6 = [cmp(6).lambda253,cmp(6).lambda273,cmp(6).lambda293,cmp(6).lambda313,cmp(6).lambda333];

lambdaliq7 = cmp(7).lambdaliq300;

%Polynomial fitting
plambda = zeros(9,2);

%plambda(1,:) = polyfit([Tliq Tgas400],[lambdaliq1 lambdagas1],1);
plambda(1,:) = polyfit(Tgas400,lambdagas1,1);
plambda(2,:) = polyfit(Tgas100,lambdagas2,1);
plambda(3,:) = polyfit(Tgas100,lambdagas3,1);
plambda(4,:) = polyfit(Tgas300,lambdagas4,1);
plambda(5,:) = polyfit(Tgas100,lambdagas5,1);
plambda(6,:) = polyfit(TLange,lambdaliq6,1);
plambda(7,:) = [0,lambdaliq7];
plambda(8,:) = [NaN,NaN];
plambda(9,:) = polyfit(Tgas200,lambdagas9,1);
end