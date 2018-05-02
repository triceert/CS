function [unt] = CAPEX_reactor(cmp,unt,str)
% INPUT: unt = unit struct
% OUTPUT: price of the reactor in US$

%% Tubes
price_tube = unt(1).price_tube;
number_of_tubes = unt(1).N_tubes;
CAPEX_reactor_priceUSdollars = price_tube*number_of_tubes;

%% Heat exchanger

%Feed
%Indexes: 1:CH4, 2:NH3
G_feed_mix = str(1).G; %total gas flow feed [mol.s-1]
y_feed = [str(1).yCH4 ; str(1).yNH3];
MW_feed = extractfield(cmp(3:4),'MW')'; %[kg.mol-1]
Tc_in = 298.15; %[K]
T = 700; %[K]
Th_in = 1600; %[K]

cp_feed = [heat_capacity(T,cmp,unt,3);heat_capacity(T,cmp,unt,4)]; %[J.mol-1.K-1]
mflow_feed = G_feed_mix*y_feed.*MW_feed; %[kg.s-1]
mflow_feed_mix = sum(mflow_feed); %[kg.s-1]

%Wall thermal diffusivity
%Assumption: 2mm wall of steel

lambdaSteel = 50.2; %[W.m-1.K-1]
deltaSteel = 0.002; %[m]
alphaSteel = lambdaSteel/deltaSteel; %[W.m-2.K-1]

%Cold feed in the heat exchanger

MW_feed_mix = sum(y_feed.*MW_feed);
cp_feed_mix = sum(y_feed.*cp_feed); % [J.mol-1.K-1]
cp_kg_feed_mix = cp_feed_mix/MW_feed_mix; % [J.kg-1.K-1]

U = alphaSteel; %[W.m-2.K-1]

%Heat flow
Q_heat = mflow_feed_mix*cp_kg_feed_mix*(T-Tc_in); %[J.s-1]

%DeltaTlm
%Assumption: at equilibrium when coming out of the heat exchanger

deltaTlm = ((Th_in - T) - (T - Tc_in))/log((Th_in - T)/(T - Tc_in)); %[K]

A_heat_exchanger = Q_heat/(U*deltaTlm); %[m2]

%Price
CAPEX_heat_exchanger_priceUSdollars = 25000*A_heat_exchanger^(0.65);%[$]

%% CAPEX in unt
unt(1).capex = CAPEX_reactor_priceUSdollars + CAPEX_heat_exchanger_priceUSdollars;
end

