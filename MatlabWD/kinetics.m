function r=kinetics(T,pCH4,pNH3,n)

r1=(7.8*1e18.*exp(-1950./T).*pCH4.*pNH3.^0.5)./...
    ((1+0.044.*exp(2390./T).*(pCH4.*pNH3.^(-0.5))).^4);
r2=(4.9*1e18.*exp(-2130./T).*pNH3)./...
    ((1+0.044.*exp(2390./T).*(pCH4.*pNH3.^(-0.5))).^4);

    switch n
        case 1
            r=r1;
        case 2
            r=r2;
    end


%antwortunits stimmen noch nicht

end