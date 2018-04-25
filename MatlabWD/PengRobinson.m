function [Z,Vm] = PengRobinson(p,T,F,cmp,unt)
% - Calculates molar volume of a compound for given pressure and
% temperature from the Peng-Robinson equation of state
%
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]       
%        F = flow [mol.s-1], vector for different components
%           cmp,unt=compund and unit struct
% OUTPUT: Vm: molar volume [m3.mol-1]
%         Z: compressibility factor


%Assign Vectors from Compound struct
%   From Nitrogen to Hydrogen Cyanide (identifier 2-6)
pc=extractfield(cmp(2:6),'pc');
Tc=extractfield(cmp(2:6),'Tc');
omega =extractfield(cmp(2:6),'omega');
R=unt(5).idgc         %[kg.m2.s-2.mol-1.K-1]



F_mix = sum(F);
z = F./F_mix; % molar fraction of each component

Tc_mix = sum(Tc.*z);
%pc_mix = sum(pc.*z);            %why we dont use this?
omega_mix = sum(omega.*z);

kappa_mix = 0.37464 + 1.54226*omega_mix - 0.26992*omega_mix^2;
Tr_mix = T./Tc_mix;
alpha_mix = (1 + kappa_mix.*(1-sqrt(Tr_mix)))^2;

a = (0.45724*R^2*Tc.^2)./pc;            
b = (0.07780*R.*Tc)./pc;
a_mix = nthroot(sum(z.*a.^(0.5)),length(a));
b_mix = sum(z.*b);

A = alpha_mix*a_mix*p/(R^2*T^2);
B = b_mix*p/(R*T);

%Polynomial solution for molar volume

a1 = p;
a2 = p*b_mix - R*T;
a3 = a_mix*alpha_mix - 3*p*b_mix^2 - 2*b_mix*R*T;
a4 = p*b_mix^3 + R*T*b_mix^2 - a_mix*b_mix*alpha_mix;

polVm = [a1 a2 a3 a4];
Vmposs = roots(polVm);

%We only want real molar volumes, and if two are there, the biggest one
Vmposs(imag(Vmposs) ~= 0) = 0;
Vm = max(Vmposs);

%Polynomial solution for compressibility factor

b1 = 1;
b2 = B-1;
b3 = A - 2*B - 3*B^2;
b4 = B^3 + B^2 - A*B;

polZ = [b1 b2 b3 b4];
Zposs = roots(polZ);

%We only want real compressibility factors
Zposs(imag(Zposs) ~= 0) = 0;
Z = Zposs(Zposs ~= 0);
end


