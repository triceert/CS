function dAdV=MBEBpfr(t,A,handles,cmp,unt,str)
 
    
    dAdV=zeros(8,1);
    %# 1 Pressure
    dAdV(1)=-handles(A(1),cmp)+3*A(3);  
    
    %#2 Nitrogen
    dAdV(2)=-A(2)+2*A(3);
    %#3 Methane
    dAdV(3)=A(1)^2-2*A(3);
    %#4 Ammonia
    dAdV(4)=0 
    %#5 Hydrogen
    dAdV(5)=0
    %#6 Hydrogen Cyanide
    dAdV(6)=0
    
    %#7 Nitrogen
    dAdV(7)=0
    %#8 Nitrogen
    dAdV(8)=0                             %

end





