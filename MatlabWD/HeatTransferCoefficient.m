function [U] = HeatTransferCoefficient(cmp,unt,p,T,F)
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        F = stream %F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%        cmp = compound struct     
%        unt = unit struct
% OUTPUT: U

R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]
MW = extractfield(cmp(1:10),'MW')';

%Wall thermal diffusivity

lambdaWall = cmp(10).lambda; %[m]
deltaWall = unt(1).deltaWall; %[m]
alphaWall = deltaWall/lambdaWall;

%In


%Out
Tout = 1600; %[K]




%U = 1/(1/alpha_in + 1/alphaWall + 1/alpha_out);

U=4.5/2.5e-3;

end

