function [plambda] = ThermalConductivity(cmp)
% - Calculates the thermal conductivity of a compound as a function of the
% temperature, calculated as a linear regression of given data
%
% INPUT: T = Temperature [K]
%        cmp = compound struct     
%
% OUTPUT: plambda 2x17 array with polynom of the compound on each line 
% Indexes: 1: Water, 2:Nitrogen,3:Methane, 4:Ammonia,5: Hydrogen, 6:Hydrogen cyanide,
% 7:Sulfuric acid, 8:Ammonium sulfate, 9:Ethane, 10:Propane, 11:Isobutane
% 12:Butane, 13:Isopentane, 14:Pentane, 15:Hexane, 16:Carbon dioxide, 17:Oxygen

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
lambdagas10 = [cmp(10).lambdagas300,cmp(10).lambdagas400,cmp(10).lambdagas500,cmp(10).lambdagas600];
lambdagas11 = [cmp(11).lambdagas300,cmp(11).lambdagas400,cmp(11).lambdagas500,cmp(11).lambdagas600];
lambdagas12 = [cmp(12).lambdagas300,cmp(12).lambdagas400,cmp(12).lambdagas500,cmp(12).lambdagas600];
lambdagas14 = [cmp(14).lambdagas400,cmp(14).lambdagas500,cmp(14).lambdagas600];
lambdagas15 = [cmp(15).lambdagas400,cmp(15).lambdagas500,cmp(15).lambdagas600];
lambdagas16 = [cmp(16).lambdagas200,cmp(16).lambdagas300,cmp(16).lambdagas400,cmp(16).lambdagas500,cmp(16).lambdagas600];
lambdagas17 = [cmp(17).lambdagas100,cmp(17).lambdagas200,cmp(17).lambdagas300,cmp(17).lambdagas400,cmp(17).lambdagas500,cmp(17s).lambdagas600];

%Other compounds
Tliq = 273.16:25:373.16;
lambdaliq1 = [cmp(1).lambdaliq273,cmp(1).lambdaliq298,cmp(1).lambdaliq323,cmp(1).lambdaliq348,cmp(1).lambdaliq373];

TLange = 253.16:20:333.16;
lambdaliq6 = [cmp(6).lambda253,cmp(6).lambda273,cmp(6).lambda293,cmp(6).lambda313,cmp(6).lambda333];

lambdaliq7 = cmp(7).lambdaliq300;
lambdaliq13 = cmp(13).lambdaliq298;

p = zeros(17,2);
p(1,:) = polyfit([Tliq Tgas400],[lambdaliq1 lambdagas1],1);
p(2,:) = polyfit(Tgas100,lambdagas2,1);
p(3,:) = polyfit(Tgas100,lambdagas3,1);
p(4,:) = polyfit(Tgas300,lambdagas4,1);
p(5,:) = polyfit(Tgas100,lambdagas5,1);
p(6,:) = polyfit(TLange,lambdaliq6,1);
p(7,:) = [0,lambdaliq7];
p(8,:) = [NaN,NaN];
p(9,:) = polyfit(Tgas200,lambdagas9,1);
p(10,:) = polyfit(Tgas300,lambdagas10,1);
p(11,:) = polyfit(Tgas300,lambdagas11,1);
p(12,:) = polyfit(Tgas300,lambdagas12,1);
p(13,:) = [0,lambdaliq13];
p(14,:) = polyfit(Tgas400,lambdagas14,1);
p(15,:) = polyfit(Tgas400,lambdagas15,1);
p(16,:) = polyfit(Tgas200,lambdagas16,1);
p(17,:) = polyfit(Tgas100,lambdagas17,1);
end