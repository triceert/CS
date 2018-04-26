function [U] = HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z)
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        F = stream %F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%        cmp = compound struct     
%        unt = unit struct
%        cp = vector (1x5) containing the different heat capacities
%        Z compressibility factor
% OUTPUT: U heat transfer coefficient

R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]
MW_in = extractfield(cmp(2:6),'MW')';
F_tot = sum(F);
D_reactor = unt(1).rad;
L = unt(1).h;

plambda = ThermalConductivity(cmp);
pmu = DynamicViscosity(cmp);

%Wall thermal diffusivity

lambdaWall = cmp(10).lambda; %[m]
deltaWall = unt(1).deltaWall; %[m]
alphaWall = deltaWall/lambdaWall;

%Inside the reactor

lambda_mix_in = MixtureThermalConductivityReactor(F,T,plambda);
mu_mix_in = MixtureDynamicViscosityReactor(F,T,pmu);
cp_mix_in = sum(F.*cp)/F_tot;
rho_mix_in = p/(R*T)*sum(F.*MW_in)/F_tot;

Pr_in = Prandtl(cp_mix_in,mu_mix_in,lambda_mix_in);

Q_in = F_tot*Z*R*T/p;
A_in = pi*D_reactor^2;
nu_mix_in = mu_mix_in/rho_mix_in;

Re_in = Reynolds(Q_in,D_reactor,A_in,nu_mix_in);

Nu_in = Nusselt_in(Re_in,Pr_in);

alpha_in = Nu_in*lambda_mix_in/D_reactor;

%Outside the reactor
%Index 1:Water, 2:CO2

y = [2;1]/3;
cp1 = heat_capacity(T,cmp,unt,1);
cp2 = heat_capacity(T,cmp,unt,9);
MW_out = [cmp(1).MW;cmp(9).MW];

lambda_mix_out = MixtureThermalConductivityNaturalGas(y,T,plambda);
mu_mix_out = MixtureDynamicViscosityNaturalGas(y,T,pmu);
cp_mix_out = y(1)*cp1 + y(2)*cp2;
%rho_mix_out = p/(R*T)*sum(y.*MW_out);

Pr_out = Prandtl(cp_mix_out,mu_mix_out,lambda_mix_out);

Re_out = 2000; %laminar

Nu_out = Nusselt_out(Re_out,Pr_out);
alpha_out = Nu_out*lambda_mix_out/L;


U = 1/(1/alpha_in + 1/alphaWall + 1/alpha_out);

%U=4.5/2.5e-3;
%U = 1/(1/alpha_in + 1/alphaWall);
end

