function [out,Ztot] = PRpartials(P,T,F,cmp,unt,n)
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
Pc=extractfield(cmp(2:6),'pc')';
Tc=extractfield(cmp(2:6),'Tc')';
w =extractfield(cmp(2:6),'omega')';
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

x=F./sum(F)


zN2=getZ(w(1),R,T,Tc(1),P*x(1),Pc(1)); 
zCH4=getZ(w(2),R,T,Tc(2),P*x(2),Pc(2)); 
zNH3=getZ(w(1),R,T,Tc(3),P*x(3),Pc(3));
zH2=getZ(w(1),R,T,Tc(4),P*x(4),Pc(4));
zHCN=getZ(w(1),R,T,Tc(5),P*x(5),Pc(5));

Ztot=(x(1).*zN2+x(2).*zCH4+x(3)*zNH3+x(4).*zH2+x(5).*zHCN);

pN2=(x(1).*zN2*P)/Ztot;
pCH4=(x(2).*zCH4*P)/Ztot;
PNH3=(x(3).*zNH3*P)/Ztot;
PH2=(x(4).*zH2*P)/Ztot;
pHCN=(x(5).*zHCN*P)/Ztot;

out=[pN2;pCH4;PNH3;PH2;pHCN];
out=out(n);

function Z=getZ(w,R,T,Tc,P,Pc)
% Reduced variables
Tr = T./Tc ;
Pc = P./Pc;

% Parameters of the EOS for a pure component
m = 0.37464 + 1.54226.*w - 0.26992.*w.^2;
alfa = (1 + m.*(1 - sqrt(Tr))).^2;
a = 0.45724.*(R.*Tc).^2./Pc.*alfa;
b = 0.0778.*R.*Tc./Pc;
A = a.*P/(R.*T).^2;
B = b.*P/(R.*T);

% Compressibility factor
Zposs = roots([1 -(1-B) (A-3.*B.^2-2.*B) -(A.*B-B.^2-B.^3)]);

 %We only want real compressibility factors
 Zposs(imag(Zposs) ~= 0) = 0;
 Z = Zposs(Zposs ~= 0);
 Z=max(Z);
 


end

end


