function [dAdV]=MBEBpfr(t,A,kinhand,parthand,cphand,Uhand,cmp,unt,str,idealreal)
%Defining odesystem
%All used handles declared in reactorptimizer.m

%% GET extra variables
%get extra variables

    H1=unt(1).deltaHHCN; %get reaction enthalpy 1
    H2=unt(1).deltaHNH; %get raction enthalpy 2
    a=unt(1).a;         %get specific surface
    Diameter=unt(1).rad*2;
    heatedwith = 0; %heated with Natural Gas: 0 / heated with H2 stream: 1
    
  
  

%F=[{'N2';'CH4';'NH3';'H2';'HCN'}]
%% INIT    
%Get PengR if idealreal set to 1 and set handle accordingly
    switch unt(1).ideal_real       
        case  1
        PRCH4=parthand(A(1),A(7),A(2:6),cmp,unt,2);
        PRNH3=parthand(A(1),A(7),A(2:6),cmp,unt,3);
        [~, Z]=parthand(A(1),A(7),A(2:6),cmp,unt,3);%get Z
        case 0            
        PRCH4=1;
        PRNH3=1;  
        Z=1;
    end
%% Assign Outputs from Handles 
    
    [r1,r2]=kinhand(A(7),A(1:6),unt,PRNH3,PRCH4,idealreal);%get rates
    cpN2=cphand(A(7),cmp,unt,2);    %get all the live cps inside reactor
    cpCH4=cphand(A(7),cmp,unt,3);
    cpNH3=cphand(A(7),cmp,unt,4);
    cpH2=cphand(A(7),cmp,unt,5);
    cpHCN=cphand(A(7),cmp,unt,6);
    
    cpH2O=cphand(A(8),cmp,unt,1);              %get all live ps outside stream
    cpCO2=cphand(A(8),cmp,unt,9) ;
    heatstreamcpmix=str(4).yH2O*cpH2O+str(4).yCO2*cpCO2;    %mixture cp of (T) for heating stream
    
    
    cprxn = [cpN2; cpCH4; cpNH3; cpH2; cpHCN];  %reaction cp
    U=Uhand(cmp,unt,A(1),A(7),A(2:6),cprxn,Z,heatedwith); %get live ht coeff
    [~,~,Re]=Uhand(cmp,unt,A(1),A(7),A(2:6),cprxn,Z,heatedwith); %get live reynolds
    [~,Q,~]=Uhand(cmp,unt,A(1),A(7),A(2:6),cprxn,Z,heatedwith);%get live vol flow 
    
    
   %other needed stuff from initial conditions
   Q_in=unt(1).Q_in;%
   rho_mix_in=unt(1).rho_mix_in;%
 
    %molar flow of heating mixture
    heatflow=str(4).G;
    
    %cocorss parameter
    cocross=unt(1).cocross;
    
%% DEFINE THE PROBLEM  ODEs  
    dAdV=zeros(9,1);
    %# 1 Pressure
    switch unt(1).pre_drop
        case 1
        dAdV(1)=-pi*(Diameter/2).^2*32*(64/Re.*(rho_mix_in.*Q_in)*Q)./(Diameter^5*pi^2);  
        case 0           
        dAdV(1)=0;
        otherwise
            disp('P drop parameter not valid')
    end
    
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
    dAdV(7)=   a *U* (A(8) - A(7)) - (r1 * H1 + r2 * H2)/...
        ( A(2) * cpN2 + A(3) * cpCH4 + A(4) * cpNH3 + A(5) * cpH2...
        + A(6) * cpHCN);
       %dAdV(7)=0;
    %#8 T heating medium
    switch unt(1).cocross
        
        case 1  
        dAdV(8)=a *U* (A(7) - A(8))./(heatstreamcpmix*heatflow);
        case 0 %crosscurrent heating
        dAdV(8)=0;      
    end
    
    %9 DQdot/Dv
    dAdV(9)=   a *U* (A(8) - A(7)) ;%



    
end





