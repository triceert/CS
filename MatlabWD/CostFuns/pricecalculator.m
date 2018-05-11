function [cmp,unt,str]=pricecalculator(cmp,unt,str)

        [unt,str]=OPEX_reactor(cmp,unt,str);
        [unt]=CAPEX_reactor(cmp,unt,str);
        [unt]=TOTEX_reactor(unt);





   %calculates overall price per k   
    tspan=10*8000*3600; %seconds operating overall ten years
    nHCN=12.86;  %moles HCN produced per second
    MWHCN=cmp(6).MW; %
    
    HCNtot=tspan*nHCN*MWHCN; %total kilogramsgrams of hcn over 10 years
    
    
    TOTEXtotal=real(unt(1).capex+unt(2).capex+unt(3).capex+unt(4).capex)*1.2+10*real(unt(1).opex+unt(2).opex+unt(3).opex+unt(4).opex);
    unt(5).toex=TOTEXtotal;
    unt(5).price=TOTEXtotal/HCNtot;
    
    %FACTOR 20% on Capex added for maintenance and personel
end
    