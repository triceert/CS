function optimizor(cmp,unt,str,idealreal,cocross)


 %% OPITMIZORRR

    %initialize
    unt(1).cocross=cocross;
    unt(1).ideal_real=idealreal;

    n=5; %optim grid size n*n

    FCH4 = linspace(1e-2,1e-4,n);    %FEEDRANGE
    ubsch=linspace(0.9,1.1,n);       %EXCESS NH3 RANGE
    nrow=linspace(1,20,n);                %NUMBER OF PFRs in row range
    hstream=linspace(1,1e-4,n);           %HEAT STREAM MOLAR FLOW
     %strin(4).G=0.1;   
    pressure=linspace(1*101325,5*101325,n)         %Pressure
    %strin(1).p=101325;
     

    %preallocate memory for different sensitivity parameters
    yieldfield=zeros(numel(n,n));
    yHCNfield=zeros(numel(n,n));
    yNH3field=zeros(numel(n,n));
    yCH4field=zeros(numel(n,n));
    pricefield=zeros(numel(n,n));

            %% FEED EXCESS PAIR
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
                     
                        %normalize results
                     yCH4field=yCH4field/(max(yCH4field(:))-min(yCH4field(:)));
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                     yHCNfield=yHCNfield/(max(yHCNfield(:))-min(yHCNfield(:)));
                     yNH3field=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));
                     yieldfield=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));
                        
                    figure
                        subplot(4,5,1)
                        contour(FCH4,ubsch,yieldfield,'ShowText','on')
                        title('Yield  resp. to CH4')

                        pbaspect([1 1 1])
                        
                        ylabel('Excess NH3')

                        subplot(4,5,2)
                        contour(FCH4,ubsch,pricefield,'ShowText','on')
                        title('Break even Price $/kg')

                        pbaspect([1 1 1])
                        

                        subplot(4,5,3)
                        contour(FCH4,ubsch,yHCNfield,'ShowText','on')
                        title('yHCN ')

                        pbaspect([1 1 1])
                        

                        subplot(4,5,4)
                        contour(FCH4,ubsch,yNH3field,'ShowText','on')
                        title('yNH3')

                        pbaspect([1 1 1])
                       

                        subplot(4,5,5)
                        contour(FCH4,ubsch,yCH4field,'ShowText','on')
                        title('yCH43')

                        pbaspect([1 1 1])
                       

                        %surf(FCH4,ubsch,yieldfield)
                        
                        
                         %% FEED Pressure PAIR
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            strin(1).p=pressure(j);  %Calculate everything           
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
                    
                         %normalize results
                     yCH4field=yCH4field/(max(yCH4field(:))-min(yCH4field(:)));
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                     yHCNfield=yHCNfield/(max(yHCNfield(:))-min(yHCNfield(:)));
                     yNH3field=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));
                     yieldfield=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));

                    subplot(4,5,6)
                    contour(FCH4,pressure,yieldfield,'ShowText','on')
                 

                    pbaspect([1 1 1])
                    
                    ylabel('Pressure')

                    subplot(4,5,7)
                    contour(FCH4,pressure,pricefield,'ShowText','on')
              

                    pbaspect([1 1 1])
               

                    subplot(4,5,8)
                    contour(FCH4,pressure,yHCNfield,'ShowText','on')
                

                    pbaspect([1 1 1])
           

                    subplot(4,5,9)
                    contour(FCH4,pressure,yNH3field,'ShowText','on')
            

                    pbaspect([1 1 1])
    

                    subplot(4,5,10)
                    contour(FCH4,pressure,yCH4field,'ShowText','on')
         

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('Pressure')
                        
                    %% FEED NROW PAIR
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
                    
                          %normalize results
                     yCH4field=yCH4field/(max(yCH4field(:))-min(yCH4field(:)));
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                     yHCNfield=yHCNfield/(max(yHCNfield(:))-min(yHCNfield(:)));
                     yNH3field=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));
                     yieldfield=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));

                    subplot(4,5,11)
                    contour(FCH4,nrow,yieldfield,'ShowText','on')
                 

                    pbaspect([1 1 1])
               
                    ylabel('PFR Segments')

                    subplot(4,5,12)
                    contour(FCH4,nrow,pricefield,'ShowText','on')
                 

                    pbaspect([1 1 1])
          

                    subplot(4,5,13)
                    contour(FCH4,nrow,yHCNfield,'ShowText','on')
                   

                    pbaspect([1 1 1])
  
                    subplot(4,5,14)
                    contour(FCH4,nrow,yNH3field,'ShowText','on')
                

                    pbaspect([1 1 1])


                    subplot(4,5,15)
                    contour(FCH4,nrow,yCH4field,'ShowText','on')
         
                    pbaspect([1 1 1])
     
                    
                    
                    %% FEED HEATMED PAIR
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            strin(4).G=hstream(j);  %Calculate everything           
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
                    
                     %normalize results
                     yCH4field=yCH4field/(max(yCH4field(:))-min(yCH4field(:)));
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                     yHCNfield=yHCNfield/(max(yHCNfield(:))-min(yHCNfield(:)));
                     yNH3field=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));
                     yieldfield=yNH3field/(max(yNH3field(:))-min(yNH3field(:)));

                    subplot(4,5,16)
                    contour(FCH4,hstream,yieldfield,'ShowText','on')
                    title('Yield with respect to CH4')

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('Molar Flow Heating Medium')

                    subplot(4,5,17)
                    contour(FCH4,hstream,pricefield,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                

                    subplot(4,5,18)
                    contour(FCH4,hstream,yHCNfield,'ShowText','on')
                  

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                

                    subplot(4,5,19)
                    contour(FCH4,hstream,yNH3field,'ShowText','on')
                 

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
     

                    subplot(4,5,20)
                    contour(FCH4,hstream,yCH4field,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
            



end