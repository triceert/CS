function [U,Q_in,Re_in] = HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z)
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        F = stream %F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%        cmp = compound struct     
%        unt = unit struct
%        cp = vector (1x5) containing the different heat capacities
%        Z compressibility factor
% OUTPUT: U heat transfer coefficient %[W.m-2.K-1]
          %Q_in flow on inside
          %Re_in reynolds on inside

R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]
MW_in = extractfield(cmp(2:6),'MW')';
F_tot = sum(F);
r_reactor = unt(1).rad;
D_reactor = 2*r_reactor;
L = unt(1).h;

%Wall thermal diffusivity

lambdaWall = cmp(10).lambda; %[W.m-1.K-1]
deltaWall = unt(1).deltaWall; %[m]
alphaWall = lambdaWall/deltaWall; %[W.m-2.K-1]

%Inside the reactor

MW_mix_in = sum(F.*MW_in)/F_tot;
lambda_mix_in = MixtureThermalConductivityReactorNew(cmp,F,T); %[W.m-1.K-1]
mu_mix_in = MixtureDynamicViscosityReactorNew(cmp,F,T); %[Pa.s]
cp_mix_in = sum(F.*cp)/F_tot; % [J.mol-1.K-1]
cp_kg_mix_in = cp_mix_in/MW_mix_in; % [J.kg-1.K-1]
rho_mix_in = Z*R*T/(p*MW_mix_in); %[m3.kg-1]

Pr_in = Prandtl(cp_kg_mix_in,mu_mix_in,lambda_mix_in);

Q_in = F_tot*Z*R*T/p; %[m3.s-1]
A_in = pi*r_reactor^2; %[m2]
nu_mix_in = mu_mix_in/rho_mix_in; %[m2.s-1]

Re_in = Reynolds(Q_in,D_reactor,A_in,nu_mix_in);

Nu_in = Nusselt_in(Re_in,Pr_in);

alpha_in = Nu_in*lambda_mix_in/D_reactor;

%Outside the reactor
%Index 1:Water, 2:CO2

y_out = [2;1]/3;
MW_out = [cmp(1).MW;cmp(9).MW];
MW_mix_out = sum(y_out.*MW_out);
cp1 = heat_capacity(T,cmp,unt,1);
cp2 = heat_capacity(T,cmp,unt,9);

lambda_mix_out = MixtureThermalConductivityNaturalGasNew(cmp,y_out,T);
mu_mix_out = MixtureDynamicViscosityNaturalGasNew(cmp,y_out,T);%/100 sollte nicht mehr nötig sein nach der Korrektur
cp_mix_out = y_out(1)*cp1 + y_out(2)*cp2; % [J.mol-1.K-1]
cp_kg_mix_out = cp_mix_out/MW_mix_out; % [J.kg-1.K-1]
%rho_mix_out = p/(R*T)*sum(y.*MW_out);

Pr_out = Prandtl(cp_kg_mix_out,mu_mix_out,lambda_mix_out);

Re_out = unt(1).Reout; %turbulent

Nu_out = Nusselt_out(Re_out,Pr_out);

V_tube = L*pi*(r_reactor + deltaWall)^2;
A_tube = L*pi*(D_reactor + 2*deltaWall);
tube_char_length = V_tube/A_tube;

alpha_out = Nu_out*lambda_mix_out/tube_char_length;




U = 1/( 1/alphaWall +1/alpha_out+1/alpha_in);



end

