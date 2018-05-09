function [opex_H2O, opex_H2SO4, opex_wastewater, opex_tot] = opex_calc(G, L, x_H2O_in, x_H2SO4_in, x_H2O_out, x_H2SO4_out, x_ammoniumsulfate_out)

time = 8000*60*60;          %time in seconds
moles_gas = time*G;         %moles needed for the molar GAS flow in one year
moles_fluid = time*L;       %moles needed for the molar FLUID flow in one year

MW_H2O = 18;    %[g/mol]
MW_H2SO4 = 98.08;   %[g/mol]
MW_ammoniumsulfat = 132.14; %

amount_water = x_H2O_in * moles_fluid * MW_H2O / 1000 / 1000;             %[ton]
amount_H2SO4 = x_H2SO4_in * moles_fluid * MW_H2SO4 / 1000 /1000;        %[ton]
amount_wastewater = (x_H2O_out*MW_H2O+x_H2SO4_out*MW_H2SO4+x_ammoniumsulfate_out*MW_ammoniumsulfat) * moles_fluid  / 1000 / 1000;     %[ton]

cost_H2O = 0.15;        %[US$/ton]
cost_H2SO4 = 100;       %[US$/ton]
cost_wastewater = 2;    %[US$/ton]



opex_H2O = cost_H2O * amount_water ;                     
opex_H2SO4 = cost_H2SO4 * amount_H2SO4 ;                 
opex_wastewater = cost_wastewater * amount_wastewater ;    
opex_tot = opex_H2O + opex_H2SO4 + opex_wastewater;             

end






