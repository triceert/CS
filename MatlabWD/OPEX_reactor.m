function [OPEX_reactor_USdollars,str] = OPEX_reactor(unt,cmp,str)
% INPUT: unt = unit struct
%        cmp = compound struct
%        str = stream struct
% OUTPUT: operating costs for the reactor [US$.s-1]
%         str = stream struct with new values for stream 2

%Constants
%Indexes: 1:CH4, 2:NH3
G_feed_mix = str(1).G; %total gas flow feed [mol.s-1]
FH2 = str(8).G*str(8).yH2; %flow of hydrogen [mol.s-1]
Qtot = unt(1).Q_tot; %total heat needed in the reactor [J.s-1]
y_feed = [str(1).yCH4 ; str(1).yNH3];
MW_feed = extractfield(cmp(3:4),'MW')'; %[kg.mol-1]
MW_CH4 = MW_feed(1);
price_natural_gas = unt(1).price_natural_gas; %[$.t-1]
deltaHc_CH4 = cmp(3).deltaHc; %[J.mol-1]
deltaHc_H2 = cmp(5).deltaHc; %[J.mol-1]

%Feeding of the reactor

price_ton = [cmp(3).price_ton ; cmp(4).price_ton]; %[$.t-1]

mflow_feed = G_feed_mix*y_feed.*MW_feed; %[kg.s-1]
price_feeding_reactor_flow = sum((mflow_feed/1000).*price_ton); %[$.s-1]

%Preheating of the stuff
%Assumption: only heated with natural gas

T = [298.15; 700]; %[K]

deltaH_CH4 = enthalpy_temperature(T(2),cmp,unt,3) - enthalpy_temperature(T(1),cmp,unt,3); %[J.mol-1]
deltaH_NH3 = enthalpy_temperature(T(2),cmp,unt,4) - enthalpy_temperature(T(1),cmp,unt,4); %[J.mol-1]
deltaH_preheating = y_feed(1)*deltaH_CH4 + y_feed(2)*deltaH_NH3; %[J.mol-1]

Q_preheat = G_feed_mix*deltaH_preheating; %[J.s-1]
FCH4_preheat = -Q_preheat/deltaHc_CH4; %[mol.s-1]
mflow_preheat_CH4 = FCH4_preheat*MW_CH4; % [kg.s-1]
price_preheating_flow_CH4 = (mflow_preheat_CH4/1000)*price_natural_gas; %[US$.s-1]

%Heating of the reactor

FCH4 = (-Qtot-FH2*deltaHc_H2)/deltaHc_CH4; %[mol.s-1]

boo = (FCH4 > 0);
switch boo
    case 0 %(FCH4 <= 0)
        price_heating_flow_CH4 = 0; %[US$.s-1]
    case 1 %(FCH4 > 0)
        mflow_CH4 = FCH4*MW_CH4; % [kg.s-1]
        price_heating_flow_CH4 = (mflow_CH4/1000)*price_natural_gas; %[US$.s-1]
end

%Total price

OPEX_reactor_USdollars = price_feeding_reactor_flow + price_preheating_flow_CH4 + price_heating_flow_CH4; %[US$.s-1]

%Values for stream 2
switch boo
    case 0 %(FCH4 <= 0)
        FCH4_tot = FCH4_preheat;
    case 1 %(FCH4 > 0)
        FCH4_tot = FCH4_preheat+FCH4;
end
str(2).G = FCH4_tot;
str(2).yCH4 = 1;

end
