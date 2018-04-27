function [Re] = Reynolds(Q,D,A,nu)
% INPUT: Q volumetric flowrate [m3.s-1]
%        D diameter [m]
%        A cross-sectional area [m2]
%        nu kinematic viscosity [m2.s-1]
% OUTPUT: Re Reynolds number
Re = Q*D/(A*nu);
end

