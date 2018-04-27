function [gamma] = vanlaar(x1, T, delta_g12, delta_g21, alpha12)
% x1 = x_LK = x_HCN
% just passing T, delta_g12, delta_g21, alpha12 since nrtl needed those 
x2 = 1-x1; % x_2 = x_HK = x_H2O
%x2 = x1;    % indices were wrong and had to relabel
%x1 = 1-x2; 
%x2=1-x1; 
A12 = 1.6418;   % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1) 
A21 = 2.7371; 
gamma(1) = exp(A12*(A21*x2/(A12*x1+A21*x2))^2); % activity coefficient of first species (HCN)
gamma(2) = exp(A21*(A12*x1/(A12*x1+A21*x2))^2); % activity coefficient of second species (H2O)
end
 