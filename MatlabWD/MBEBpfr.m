function dAdV=MBEBpfr(t,A,kinhand,parthand,cphand,cmp,unt,str,idealreal)
%Defining odesystem
%All used handles declared in reactorptimizer.m

%% Provisorisch 
    H1=unt(1).deltaHHCN;
    H2=unt(1).deltaHNH; %does not load
    %H2=91880;

  
    a=unt(1).a;
    U=4.5/2.5e-3;
    cp=71000;

%F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%% INIT    
%Get PengR if idealreal set to 1 and set handle accordingly
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
    cpN2=cphand(A(7),cmp,unt,2)
    cpCH4=cphand(A(7),cmp,unt,3)
    cpNH3=cphand(A(7),cmp,unt,4)
    cpH2=cphand(A(7),cmp,unt,5)
    cpHCN=cphand(A(7),cmp,unt,6)
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
    dAdV(7)=   a *(U* (A(8) - A(7)) - (r1 * H1 + r2 * H2))/...
        ( A(2) * cpN2 + A(3) * cpCH4 + A(4) * cpNH3 + A(5) * cpH2...
        + A(6) * cpHCN);
       %dAdV(7)=0;
    %#8 T heating medium
    dAdV(8)=0;                             %


    
    
    
end





