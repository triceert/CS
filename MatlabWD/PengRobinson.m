function [Vmposs] = PengRobinson(p,T,pc,Tc,omega)
% - Calculates molar volume of a compound for given pressure and
% temperature from the Peng-Robinson equation of state
%
% INPUT: p = pressure  [Pa]
%        T = Temperature [K]
%        pc = critical pressure [Pa]
%        Tc = critical temperature [K]
%        omega = acentric factor

R = 8.314; % [kg.m2.s-2.mol-1.K-1]
kappa = 0.37464 + 1.54226*omega - 0.26992*omega^2;
Tr = T/Tc;
alpha = (1 + kappa*(1-sqrt(Tr)))^2;
a = (0.45724*R^2*Tc^2)/pc;
b = (0.07780*R*Tc)/pc;

%Polynomial solution

a1 = p;
a2 = p*b - R*T;
a3 = a*alpha - 3*p*b^2 - 2*b*R*T;
a4 = p*b^3 + R*T*b^2 - a*b*alpha;

delta0 = a2^2 - 3*a1*a3;
delta1 = 2*a2^3 - 9*a1*a2*a3 + 27*a1^2*a4;
C(1) = nthroot((delta1 + sqrt(delta1^2 - 4*delta0^3))/2,3);
C(2) = C(1)*(-0.5+0.5*i*sqrt(3));
C(3) = C(1)*(-0.5-0.5*i*sqrt(3));

Vmposs(1) = -1/(3*a)*(b+C(1)+delta0/C(1));
Vmposs(2) = -1/(3*a)*(b+C(2)+delta0/C(2));
Vmposs(3) = -1/(3*a)*(b+C(3)+delta0/C(3));

%We only want real molar volumes
end


