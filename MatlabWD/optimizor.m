function optimizor(cmp,unt,str,idealreal,cocross)


 %% OPITMIZORRR

    %initialize
    unt(1).cocross=cocross;
    unt(1).ideal_real=idealreal;

    n=5; %optim grid size n*n

    FCH4 = linspace(0.1,0.0001,n);    %FEEDRANGE
    ubsch=linspace(0.9,1.1,n);       %EXCESS NH3 RANGE
    nrow=linspace(1,10,n);                %NUMBER OF PFRs in row range

    %preallocate memory
    yieldfield=zeros(numel(n,n));
    %ntubesfield=zeros(numel(n,n));
    yHCNfield=zeros(numel(n,n));
    yNH3field=zeros(numel(n,n));
    yCH4field=zeros(numel(n,n));
    pricefield=zeros(numel(n,n));

            %FEED EXCESS PAIR
                        for i=1:n

                            for j=1:n
                                str(1).FCH4=FCH4(i);
                                str(1).ubsch=ubsch(j);  %Calculate everything           
                                        [cmp,unt,str]=reactorcalculator(cmp,unt,str,0);   %plotparameter 0
                                        [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
                                        [cmp,unt,str] = hcnabsorption2(cmp,unt,str);       
                                        [cmp,unt,str]=hcn_distillation(cmp,unt,str); 
                                        [unt,str]=OPEX_reactor(cmp,unt,str);
                                        [unt]=CAPEX_reactor(unt);
                                        [unt]=TOTEX_reactor(unt);

                                pricefield(i,j)=pricecalculator(unt,cmp);
                                yieldfield(i,j)=real(unt(1).yield);                                                        
                                yHCNfield(i,j)=real(str(5).yHCN);
                                yNH3field(i,j)=real(str(5).yNH3);
                                yCH4field(i,j)=real(str(5).yCH4);

                            end

                        end




                        figure
                        subplot(3,5,1)
                        contour(FCH4,ubsch,yieldfield,'ShowText','on')
                        title('Yield with respect to CH4')

                        pbaspect([1 1 1])
                        xlabel('Feed CH4')
                        ylabel('Excess NH3')

                        subplot(3,5,2)
                        contour(FCH4,ubsch,pricefield,'ShowText','on')
                        title('Break even Price $/kg')

                        pbaspect([1 1 1])
                        xlabel('Feed CH4')
                        ylabel('Excess NH3')

                        subplot(3,5,3)
                        contour(FCH4,ubsch,yHCNfield,'ShowText','on')
                        title('yHCN ')

                        pbaspect([1 1 1])
                        xlabel('Feed CH4')
                        ylabel('Excess NH3')

                        subplot(3,5,4)
                        contour(FCH4,ubsch,yNH3field,'ShowText','on')
                        title('yNH3')

                        pbaspect([1 1 1])
                        xlabel('Feed CH4')
                        ylabel('Excess NH3')

                        subplot(3,5,5)
                        contour(FCH4,ubsch,yCH4field,'ShowText','on')
                        title('yCH43')

                        pbaspect([1 1 1])
                        xlabel('Feed CH4')
                        ylabel('Excess NH3')

                        %surf(FCH4,ubsch,yieldfield)
                        
                        
                    %FEED NROW PAIR
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            unt(1).nrow=nrow(j);  %Calculate everything           
                                    [cmp,unt,str]=reactorcalculator(cmp,unt,str,0);   %plotparameter 0
                                    [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
                                    [cmp,unt,str] = hcnabsorption2(cmp,unt,str);       
                                    [cmp,unt,str]=hcn_distillation(cmp,unt,str); 
                                    [unt,str]=OPEX_reactor(cmp,unt,str);
                                    [unt]=CAPEX_reactor(unt);
                                    [unt]=TOTEX_reactor(unt);

                            pricefield(i,j)=pricecalculator(unt,cmp);
                            yieldfield(i,j)=real(unt(1).yield);                                                        
                            yHCNfield(i,j)=real(str(5).yHCN);
                            yNH3field(i,j)=real(str(5).yNH3);
                            yCH4field(i,j)=real(str(5).yCH4);

                        end

                    end



                    
                    
                    subplot(3,5,6)
                    contour(FCH4,nrow,yieldfield,'ShowText','on')
                    title('Yield with respect to CH4')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('PFR Segments')

                    subplot(3,5,7)
                    contour(FCH4,nrow,pricefield,'ShowText','on')
                    title('Break even Price $/kg')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('PFR Segments')

                    subplot(3,5,8)
                    contour(FCH4,nrow,yHCNfield,'ShowText','on')
                    title('yHCN ')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('PFR Segments')

                    subplot(3,5,9)
                    contour(FCH4,nrow,yNH3field,'ShowText','on')
                    title('yNH3')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('PFR Segments')

                    subplot(3,5,10)
                    contour(FCH4,nrow,yCH4field,'ShowText','on')
                    title('yCH43')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('PFR Segments')



end