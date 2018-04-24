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

F_tot = sum(F);
z = F./F_tot; % molar fraction of each component

kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
Tr = T./Tc;
alpha = (1 + kappa.*(1-sqrt(Tr))).^2;

a = (0.45724*R^2*Tc.^2)./pc;
b = (0.07780*R.*Tc)./pc;
a_tot = 0;
b_tot = 0;
for j = 1:length(a)
    a_tot = a_tot + z(j)*sqrt(a(j));
    b_tot = b_tot + z(j)*b(j);
end
a_tot = nthroot(a_tot,length(a));

A = alpha*a_tot*p/(R^2*T^2);
B = b_tot*p/(R*T);

%Polynomial solution for molar volume

a1 = p;
a2 = p*b_tot - R*T;
a3 = a_tot*alpha - 3*p*b_tot^2 - 2*b_tot*R*T;
a4 = p*b_tot^3 + R*T*b_tot^2 - a_tot*b_tot*alpha;

delta0Vm = a2^2 - 3*a1*a3;
delta1Vm = 2*a2^3 - 9*a1*a2*a3 + 27*a1^2*a4;
CVm(1) = nthroot((delta1Vm + sqrt(delta0Vm^3))/2,3);
CVm(2) = CVm(1)*(-0.5+0.5*i*sqrt(3));
CVm(3) = CVm(1)*(-0.5-0.5*i*sqrt(3));

for k = 1:3
    Vmposs(k) = -1/(3*a_tot)*(b_tot+CVm(k)+delta0Vm/CVm(k));
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
    Zposs(k) = -1/(3*a_tot)*(b_tot+CZ(k)+delta0Z/CZ(k));
end

%We only want real compressibility factors
Zposs(imag(Vmposs) ~= 0) = 0;
Z = Zposs(Zposs ~= 0));
end


