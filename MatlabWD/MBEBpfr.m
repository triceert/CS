function dAdV=MBEBpfr(t,A,kinhand,parthand,cmp,unt,str)
 
%kinhand = @(T,F), -- nth reaction
%parthand  = @(p,T,F,cmp,unt,n) -- partial pressure of nth comp
%F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
    

    [r1,r2]=kinhand(A(7),A(1:6));
  


    dAdV=zeros(8,1);
    %# 1 Pressure
    dAdV(1)=0;  
    
    %#2 Nitrogen
    dAdV(2)=r2;  
    %#3 Methane
    dAdV(3)=-r1;
    %#4 Ammonia
    dAdV(4)=(-r1-2*r2);
    %#5 Hydrogen
    dAdV(5)=(3*r1+3*r2);
    %#6 Hydrogen Cyanide
    dAdV(6)=r1;
    
    %#7 Trxn
    dAdV(7)=0;
    %#8 T heating medium
    dAdV(8)=0;                             %

    %P=A(2:6)
    
    
    
end





