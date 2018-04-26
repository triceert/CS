function dAdV=MBEBpfr(t,A,kinhand,parthand,cmp,unt,str,idealreal)
%% Provisorisch 
    H1=unt(1).deltaHHCN;
    %H2=unt(2).deltaHNH3 does not load
    H2=91880;

  
    a=unt(1).a;
    U=4.5/2.5e-3;
    cp=71000;

%F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%% INIT    
%Get PR if idealreal set to 1
    switch idealreal       
        case  1
        PRCH4=parthand(A(1),A(7),A(2:6),cmp,unt,2);
        PRNH3=parthand(A(1),A(7),A(2:6),cmp,unt,3);        
        case 0            
        PRCH4=1;
        PRNH3=1;            
    end
%% Assign Outputs from Handles    
    [r1,r2]=kinhand(A(7),A(1:6),unt,PRNH3,PRCH4,idealreal);

%% DEFINE THE PROBLEM  ODEs  
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
    %dAdV(7)=   a *(U* (A(8) - A(7)) - (r1 * H1 + r2 * H2))/...
       % ( A(2) * cp + A(3) * cp + A(4) * cp + A(5) * cp...
       % + A(6) * cp);
       dAdV(7)=0;
    %#8 T heating medium
    dAdV(8)=0;                             %


    
    
    
end





