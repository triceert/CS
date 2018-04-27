function [r1,r2]=kinetics(T,F,unt,PRNH3,PRCH4,idealreal)

    %idealreal 0 ideal 1 real(Peng robinson)


    
    Factor=(unt(1).a*1e4)/6.022e23; 
    %Correction Factor with specific surface
    %Avogadro und Flächenkorrektur für 
    

%% IDEAL REAL GAS SWITCH

    switch idealreal
    
        case 0
        pCH4=(F(3)/sum(F(2:6))).*F(1);
        pNH3=(F(4)/sum(F(2:6))).*F(1);
    
        case 1
        pCH4=PRCH4*(F(3)/sum(F(2:6))).*F(1);
        pNH3=PRNH3*(F(4)/sum(F(2:6))).*F(1);
        
        otherwise
            error('Wrong idealreal-parameter')
            
    end
    

    %ptto=0.00750062; %ptotorr factor
    pCH4=0.00750062*pCH4;   %AND GET THE FUCKING IMPERIAL SHIT TOGETHER  
    pNH3=0.00750062*pNH3; 

%%
%Get Kinetic constants

r1=Factor.*(7.8*1e18*exp(-1950./T)*pCH4*pNH3.^0.5)/...
    ((1+0.044*exp(2390./T)*pCH4.*pNH3^(-0.5))^4);
r2=Factor .*(4.9*1e18*exp(-2130./T)*pNH3)/...
    ((1+0.044*exp(2390./T)*pCH4*pNH3^(-0.5))^3);

    



end