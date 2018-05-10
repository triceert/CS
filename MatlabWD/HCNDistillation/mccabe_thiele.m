function = mccabe_thiele()
x_plot = 0:0.01:1;
y_plot = x_plot; 
delta_g12 = 1298.9610; % data from J. Gmehling, U. Onken, W. Arlt, Vapor-Liquid Equilibrium Data Collection, Aqueous-Organic Systems (Supplement 1)
delta_g21 = 539.9577; 
alpha12 = 0.3836; % this has nothing to do with the relative volatility alpha
nrtl(x1, , delta_g12, delta_g21, alpha12)





end