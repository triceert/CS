function [U] = HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z)
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        F = stream %F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%        cmp = compound struct     
%        unt = unit struct
%        cp = vector (1x5) containing the different heat capacities
%        Z compressibility factor
% OUTPUT: U

R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]
MW = extractfield(cmp(1:10),'MW')';
F_tot = sum(F);
D_reactor = unt(1).rad;

%Wall thermal diffusivity

lambdaWall = cmp(10).lambda; %[m]
deltaWall = unt(1).deltaWall; %[m]
alphaWall = deltaWall/lambdaWall;

%In

plambda_in = ThermalConductivity(cmp);
pmu_in = DynamicViscosity(cmp);

F_in = [0 F 0 0 0];
lambda_mix_in = MixtureThermalConductivity(F_in,T,plambda_in);
mu_mix_in = MixtureDynamicViscosity(F_in,T,pmu_in);
cp_mix_in = sum(F.*cp)/F_tot;
rho_mix_in = p/(R*T)*(F_in.*MW);

Pr_in = Prandtl(cp_mix_in,mu_mix_in,lambda_mix_in);

Q_in = F_tot*Z*R*T/p;
A_in = pi*D_reactor^2;
nu_mix_in = mu_mix_in/rho_mix_in;

Re_in = Reynolds(Q_in,D_reactor,A_in,nu_mix_in);

Nu_in = Nusselt_in(Re_in,Pr_in);

alpha_in = Nu_in*lambda_mix_in/D_reactor;

%Out

%Tout = 1600; %[K]




%U = 1/(1/alpha_in + 1/alphaWall + 1/alpha_out);

%U=4.5/2.5e-3;
U = 1/(1/alpha_in + 1/alphaWall);
end

