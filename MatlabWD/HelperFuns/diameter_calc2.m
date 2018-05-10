function diameter = diameter_calc2(velocity,G, y_H2_in, y_N2_in, y_HCN_in, y_CH4_in, y_NH3_in);

M_H2 = 2;      % Molar masses of the compounds in [g/mol]
M_N2 = 28;
M_HCN = 27;
M_CH4 = 16;
M_NH3 = 17;

velocity_ft = 3.2808* velocity;   % velocity in [ft/s]
vapor_dens = (1/velocity_ft)^2;    % vapor density in [lb/ft^3]
vapor_dens_si = vapor_dens*453.59237/0.028317;      % vapor density in [g/m^3]
m_flow = G * (y_H2_in * M_H2 + y_N2_in * M_N2 + y_HCN_in * M_HCN + y_CH4_in * M_CH4 + y_NH3_in * M_NH3);      % mass flow in [g/s]
a_cross = m_flow/(velocity*vapor_dens_si);      % Cross sectional area in [m^2]

diameter = sqrt((4 * a_cross)/pi);        % Diameter in m

end