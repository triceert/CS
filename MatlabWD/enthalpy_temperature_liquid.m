function [HT_L] = enthalpy_temperature_liquid(T,cmp,unt)
enthalpy_vapor_at_boiling_point = zeros(6,1); 
enthalpy_liquid_at_boiling_point = zeros(6,1); 
HT_L = zeros(6,1); 
for i=1:6
    HT(i) = enthalpy_temperature(cmp(i).bp,cmp,unt, i);
    %enthalpy_vapor_at_boiling_point(i) = HT(i); % need to extract only that one component for which we have just used the boiling point
    clear HT;
    enthalpy_liquid_at_boiling_point(i) = enthalpy_vapor_at_boiling_point(i) - cmp(i).Hv; 
    HT_L(i) = enthalpy_liquid_at_boiling_point(i) -cmp(i).cpl.*(cmp(i).bp-T); 
end

