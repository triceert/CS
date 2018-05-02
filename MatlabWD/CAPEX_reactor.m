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
lambda_feed_mix = MixtureThermalConductivityHeatExchanger(cmp,y_feed,T); %[W.m-1.K-1]
mu_feed_mix = MixtureDynamicViscosityReactorHeatExchanger(cmp,y_feed,T); %[Pa.s]
cp_feed_mix = sum(y_feed.*cp_feed); % [J.mol-1.K-1]
cp_kg_feed_mix = cp_feed_mix/MW_feed_mix; % [J.kg-1.K-1]

Pr_feed = Prandtl(cp_kg_feed_mix,mu_feed_mix,lambda_feed_mix);
Re_feed = 10000; %turbulent
Nu_feed = Nusselt_in(Re_feed,Pr_feed);

%Natural gas heating
%Index 1:Water, 2:CO2

Th_in = 1600; %[K]

y_NG = [2;1]/3;
MW_NG = [cmp(1).MW;cmp(9).MW];
MW_mix_NG = sum(y_NG.*MW_NG);
cp1 = heat_capacity(T,cmp,unt,1);
cp2 = heat_capacity(T,cmp,unt,9);

lambda_mix_NG = MixtureThermalConductivityNaturalGasNew(cmp,y_NG,T); %[W.m-1.K-1]
mu_mix_NG = MixtureDynamicViscosityNaturalGasNew(cmp,y_NG,T); %[Pa.s]
cp_mix_NG = y_NG(1)*cp1 + y_NG(2)*cp2; % [J.mol-1.K-1]
cp_kg_mix_NG = cp_mix_NG/MW_mix_NG; % [J.kg-1.K-1]

Pr_NG = Prandtl(cp_kg_mix_NG,mu_mix_NG,lambda_mix_NG);
Re_NG = unt(1).Reout; %turbulent
Nu_NG = Nusselt_out(Re_NG,Pr_NG);

%Heat flow
Q_heat = mflow_feed_mix*cp_kg_feed_mix*(T-Tc_in); %[J.s-1]

%DeltaTlm
%Assumption: at equilibrium when coming out of the heat exchanger

deltaTlm = ((Th_in - T) - (T - Tc_in))/log((Th_in - T)/(T - Tc_in)); %[K]

%Finding the diameter D of the heat exchanger

a = pi/4;
b = -(Q_heat/deltaTlm)*(1/(Nu_feed*lambda_feed_mix)+1/(Nu_NG*lambda_mix_NG));
c = -(Q_heat/deltaTlm)/alphaSteel;

polD = [a b c];
Dposs = roots(polD);
%We only want real molar volumes bigger than 0
Dposs(imag(Dposs) ~= 0) = 0;
D_heat_exchanger = max(Dposs);

A_heat_exchanger = pi*(D_heat_exchanger/2)^2; %[m2]

%Price
CAPEX_heat_exchanger_priceUSdollars = 25000*A_heat_exchanger^(0.65)%[$]

%% CAPEX in unt
unt(1).capex = CAPEX_reactor_priceUSdollars + CAPEX_heat_exchanger_priceUSdollars;
end

