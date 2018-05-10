function [cmp,unt,str] = hcnabsorber(cmp,unt,str)


    switch unt(1).ideal_real
        case 0 %ideal case
            [cmp, unt, str] = hcnideal(cmp, unt, str);
            
        case 1 %real case
            [cmp, unt, str] = hcnnonideal(cmp, unt, str);
            
    end


end