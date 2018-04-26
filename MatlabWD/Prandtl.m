function [Pr] = Prandtl(cp,mu,lambda)
% INPUT: cp heat capacity [J.mol-1.K-1]
%        mu dynamic viscosity [Pa.s]
%        lambda thermal conductivity [W.m-1.K-1]
% OUTPUT: Pr Prandtl number
Pr = cp*mu/lambda;
end

