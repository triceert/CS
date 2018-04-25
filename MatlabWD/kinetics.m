function [r1,r2]=kinetics(T,F)



ptto=0.00750062; %ptotorr factor


    specvol=2.667; %cm2/cm3
    Factor=specvol/6.022e23;
    Factor=Factor.*1e4 ;%to get to m2




%for idg
pCH4=(F(3)/sum(F(2:6))).*F(1);
pNH3=(F(4)/sum(F(2:6))).*F(1);




pCH4=ptto*pCH4;
pNH3=ptto*pNH3;



r1=Factor.*(7.8*1e18*exp(-1950./T)*pCH4*pNH3.^0.5)/...
    ((1+0.044*exp(2390./T)*pCH4.*pNH3^(-0.5))^4);
r2=Factor .*(4.9*1e18*exp(-2130./T)*pNH3)/...
    ((1+0.044*exp(2390./T)*pCH4*pNH3^(-0.5))^3);

    



end