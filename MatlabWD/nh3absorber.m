function [cmp,unt,str] = nh3absorber(cmp,unt,str)


    switch unt(1).ideal_real
        case 0 %ideal case
            [cmp, unt, str] = NH3_absorber_ideal(cmp, unt, str);
            
        case 1 %real case
            [cmp, unt, str] = NH3_absorber_nonideal4(cmp, unt, str);
            
    end


end