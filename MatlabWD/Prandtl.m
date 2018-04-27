function [Pr] = Prandtl(cp,mu,lambda)
% INPUT: cp heat capacity [J.kg-1.K-1] !!!!!!!!!!!!!!!!!! kg and not mole,
% as calculated in the heat capacity
%        mu dynamic viscosity [Pa.s]
%        lambda thermal conductivity [W.m-1.K-1]
% OUTPUT: Pr Prandtl number
Pr = cp*mu/lambda;
end

