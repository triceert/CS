function [Vm,Z] = PengRobinson(p,T,pc,Tc,omega,F)
% - Calculates molar volume of a compound for given pressure and
% temperature from the Peng-Robinson equation of state
%
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        pc = critical pressure [Pa], vector for different components
%        Tc = critical temperature [K], vector for different components
%        omega = acentric factor, vector for different components
%        F = flow [mol.s-1], vector for different components
% OUTPUT: Vm: molar volume [m3.mol-1]
%         Z: compressibility factor

R = 8.314; % [kg.m2.s-2.mol-1.K-1]

F_mix = sum(F);
z = F./F_mix; % molar fraction of each component

Tc_mix = sum(Tc.*z);
pc_mix = sum(pc.*z);
omega_mix = sum(omega.*z);

kappa_mix = 0.37464 + 1.54226*omega_mix - 0.26992*omega_mix^2;
Tr_mix = T/Tc_mix;
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

delta0Vm = a2^2 - 3*a1*a3;
delta1Vm = 2*a2^3 - 9*a1*a2*a3 + 27*a1^2*a4;
CVm(1) = nthroot((delta1Vm + sqrt(delta0Vm^3))/2,3);
CVm(2) = CVm(1)*(-0.5+0.5*i*sqrt(3));
CVm(3) = CVm(1)*(-0.5-0.5*i*sqrt(3));

for k = 1:3
    Vmposs(k) = -1/(3*a_mix)*(b_mix+CVm(k)+delta0Vm/CVm(k));
end

%We only want real molar volumes, and if two are there, the biggest one
Vmposs(imag(Vmposs) ~= 0) = 0;
Vm = max(Vmposs);

%Polynomial solution for compressibility factor

b1 = 1;
b2 = B-1;
b3 = A - 2*B - 3*B^2;
b4 = B^3 + B^2 - A*B;

delta0Z = b2^2 - 3*b1*b3;
delta1Z = 2*b2^3 - 9*b1*b2*b3 + 27*b1^2*b4;
CZ(1) = nthroot((delta1Z + sqrt(delta1Z^2 - 4*delta0Z^3))/2,3);
CZ(2) = CZ(1)*(-0.5+0.5*i*sqrt(3));
CZ(3) = CZ(1)*(-0.5-0.5*i*sqrt(3));

for k = 1:3
    Zposs(k) = -1/(3*a_mix)*(b_mix+CZ(k)+delta0Z/CZ(k));
end

%We only want real compressibility factors
Zposs(imag(Vmposs) ~= 0) = 0;
Z = Zposs(Zposs ~= 0));
end


