function [out,ZTot] = PRpartials(P,T,F,cmp,unt,n)
% - Calculates partial pressure of nth component of F of a compound for given pressure and
% temperature from the Peng-Robinson equation of state
%
% Second output, Z of PFR reaction mixture

% INPUT: p = pressure  [Pa]
%        T = Temperature [K]       
%        F = flow [mol.s-1], vector for different components
%           cmp,unt=compund and unit struct
%        n=n th vcompenent of F vector; n is useless if one is just
%        interessted in second output Z
% OUTPUT: pi=partial pressure
%         


%Assign Vectors from Compound struct
%   From Nitrogen to Hydrogen Cyanide (identifier 2-6)
% Pc=extractfield(cmp(2:6),'pc')';
% Tc=extractfield(cmp(2:6),'Tc')';
% w =extractfield(cmp(2:6),'omega')';
% R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]

pc=extractfield(cmp(2:6),'pc')';
Tc=extractfield(cmp(2:6),'Tc')';
omega =extractfield(cmp(2:6),'omega')';
R=unt(5).idgc;         %[kg.m2.s-2.mol-1.K-1]







% Symbols:
% P = final pressure of the gas mixture.
% pza, pzb, pzc, ... pzn = corrected partial pressure of the gas.
% a, b, c, ...., n = concentration of the specific gas (%Gas/100).
% za, zb, zc, .., zn = compressibility factors of the different components at pressures a*P, b*P, c*P, ...n*P. Ztot = compressibility factor for the gas mixture, and
% Ztot =(a*za)+(b*zb)+(c*zc)+....+(n*zn)
% pza =(a*za *P)/Ztot
% pzb =(b*zb *P)/Ztot
% pzc =(c*zc *P)/Ztot

%[{'N2';'CH4';'NH3';'H2';'HCN'}]

x=F./sum(F);


 % This function calculates the parameter for the peng robinson equation
 %1:NH3
 %2:CH4
 %3:HCN
 %4:H2
 %5:N2
 % T: temperature in K
%577 % x: mole fractions (dimensionless)
%579 R = 0.062363577;
% 580 P=P0;
% 581 Tc = [132.25, -82.59, 183.5, -240.01, -146.96] + 273.15;
% 582 pc = [113.3, 45.99, 53.9, 12.96, 33.96] .* 750.06;
%(mˆ3*torr)/(mol K)
%ideal pressure in torr %critical temperatures in K
%critical pressures in torr
%reduced temperatures (dimensionless)
 Tr=T./Tc;
 s = [0.37464 ,1.54226 , -0.26992];
%omega = [0.25, 0.011, 0.4095, -0.216, 0.039];


omegamix = x(1) .* omega(1) + x(2) .* omega(2) + x(3) .* omega(3) + ...
 x(4) .* omega(4) + x(5) .* omega(5);
 S = s(1) + s(2) .* omegamix + s(3) .* omegamix.^2;
 k=(1+S.*(1-sqrt(Tr))).^2;
 E=S.*sqrt(Tr./k);
 a=(0.457235.*(R).*(Tc.^2).*k)./pc;
 b = (0.07780 .* R .* Tc) ./ (pc);
 amix = (x(1) .* sqrt(a(1)) + x(2) .* sqrt(a(2)) + ...
 x(3) .* sqrt(a(3)) + x(4) .* sqrt(a(4)) + ...
 x(5) .* sqrt(a(5))).^2;
 bmix= x(1).*b(1)+x(2).*b(2)+x(3).*b(3)+...
 x(4) .* b(4) + x(5) .* b(5);
 Amix=(amix.*P)./((R.*T).^2);
 Bmix=(bmix.*P)./(R.*T);
  
coef = [1, (-1 + Bmix), (Amix - 2 * Bmix - 3 * Bmix^2), (-Amix * Bmix + Bmix^2 + Bmix^3)];
Z = roots(coef);
Z = max(Z);


ZTot=Z;

% ^
% 
% zN2=getZ(w(1),R,T,Tc(1),P*x(1),Pc(1)); 
% zCH4=getZ(w(2),R,T,Tc(2),P*x(2),Pc(2)); 
% zNH3=getZ(w(3),R,T,Tc(3),P*x(3),Pc(3));
% zH2=getZ(w(4),R,T,Tc(4),P*x(4),Pc(4));
% zHCN=getZ(w(5),R,T,Tc(5),P*x(5),Pc(5));
% 
% Ztot=(x(1).*zN2+x(2).*zCH4+x(3)*zNH3+x(4).*zH2+x(5).*zHCN)
% 
% pN2=(x(1).*zN2*P)/Ztot;
% pCH4=(x(2).*zCH4*P)/Ztot;
% PNH3=(x(3).*zNH3*P)/Ztot;
% PH2=(x(4).*zH2*P)/Ztot;
% pHCN=(x(5).*zHCN*P)/Ztot;
% 
% out=[pN2;pCH4;PNH3;PH2;pHCN];
% out=out(n);
% 
% function Z=getZ(w,R,T,Tc,P,Pc)
% % Reduced variables
% Tr = T./Tc ;
% Pc = P./Pc;
% 
% % Parameters of the EOS for a pure component
% m = 0.37464 + 1.54226.*w - 0.26992.*w.^2;
% alfa = (1 + m.*(1 - sqrt(Tr))).^2;
% a = 0.45724.*(R.*Tc).^2./Pc.*alfa;
% b = 0.0778.*R.*Tc./Pc;
% A = a.*alfa.*P/(R.^2*T^2);
% B = b.*P/(R.*T);
% 
% % Compressibility factor
% Zposs = roots([1 -(1-B) (A-3.*B.^2-2.*B) -(A.*B-B.^2-B.^3)]);
% 
%  %We only want real compressibility factors
%  Zposs(imag(Zposs) ~= 0) = 0;
%  Z = Zposs(Zposs ~= 0);
%  Z=Z(end);
%  
% 
% 
% end
Videal=R.*T/P;
Vreal=Z.*Videal;

Preal=R.*T/Vreal;
out=Preal.*x(n);
out=P.*x(n);
end




