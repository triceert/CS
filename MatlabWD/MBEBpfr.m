function dAdV=MBEBpfr(t,A,kinhand,prhand,cmp,unt,str)
 
%pch4handle
%pnh3handle


    
    dAdV=zeros(8,1);
    %# 1 Pressure
    dAdV(1)=kinhand(200,101354,101354,1);  
    
    %#2 Nitrogen
    dAdV(2)=kinhand(200,101354,101354,2);  
    %#3 Methane
    dAdV(3)=0;
    %#4 Ammonia
    dAdV(4)=0; 
    %#5 Hydrogen
    dAdV(5)=0;
    %#6 Hydrogen Cyanide
    dAdV(6)=0;
    
    %#7 Nitrogen
    dAdV(7)=0;
    %#8 Nitrogen
    dAdV(8)=0;                             %

end





