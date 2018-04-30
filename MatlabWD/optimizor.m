function optimizor(cmp,unt,str)


 %% OPITMIZORRR

    %initialize
    disp('Optimization started, dependent on the performance of your computer this can take up to 5 minutes. Its a good time to grab a coffe now.')

    n=4; %optim grid size n*n
    
    %get persistent variables
    persFCH4=str(1).FCH4;
    persRatio=str(1).ubsch;
    persPressure=str(1).p;
    persnrow=unt(1).nrow;

    FCH4 = linspace(5e-2,1e-4,n);    %FEEDRANGE
    ubsch=linspace(0.9,1.1,n);       %EXCESS NH3 RANGE
    nrow=linspace(1,10,n);                %NUMBER OF PFRs in row range
    hstream=linspace(1e-1,1e-4,n);           %HEAT STREAM MOLAR FLOW
    pressure=linspace(1*101325,3*101325,n);         %Pressure
    temperature=linspace(600,800,n);
     

    %preallocate memory for different sensitivity parameters
    yieldfield=zeros(numel(n,n));
    yHCNfield=zeros(numel(n,n));
    yNH3field=zeros(numel(n,n));
    yH2field=zeros(numel(n,n));
    pricefield=zeros(numel(n,n));
    
    
    
%     
%     function price=optimfunc(x,cmp,unt,str)
%     
%         
%           str(1).FCH4=x(1);
%           str(1).ubsch=x(2);
%           str(1).p=x(3);
%           str(1).T=x(4);
%           unt(1).nrow=x(5);       
%           str(4).G=x(6);
%          
%           
%        
%                             [cmp,unt,str]=reactorcalculator(cmp,unt,str,0);   %plotparameter 0
%                                         [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
%                                         [cmp,unt,str] = hcnabsorption2(cmp,unt,str);       
%                                         [cmp,unt,str]=hcn_distillation(cmp,unt,str); 
%                                         [unt,str]=OPEX_reactor(cmp,unt,str);
%                                         [unt]=CAPEX_reactor(unt);
%                                         [unt]=TOTEX_reactor(unt);
% 
%                                 price=pricecalculator(unt,cmp);
%                                 
%     end
% 
%      optimhandle=@(x)optimfunc(x,cmp,unt,str)
%      options.Algorithm = 'levenberg-marquardt';
%      
%       lb = [-Inf,-Inf,-Inf,-Inf,-Inf,-Inf];
%       ub = [Inf, Inf, Inf, Inf, Inf, Inf];
%       lb=[];
%       ub=[];
%       x0 = [0.018,1.05,101325,700,1,0.01];
%       x = lsqnonlin(optimhandle,x0,lb,ub,options)
% 
%         optimfeed=x(1)
%         optiexcess=x(2)
%         optipressure=x(3)
%         optitemperature=x(4)
%         optirow=x(5)
%         optiflow=x(6)
      
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
                                yH2field(i,j)=real(str(5).yH2);

                            end

                        end
                     
                     %normalize results
                     
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                    
                     
                     %find min  or max
                     
                        
                    figure
                        subplot(5,5,1)
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yieldfield,'ShowText','on')
                        title('Yield  resp. to CH4')

                        pbaspect([1 1 1])
                        
                        ylabel('Excess NH3')

                        subplot(5,5,2)
                                                hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,pricefield,'ShowText','on')
                        title('Break even Price $/kg')

                        pbaspect([1 1 1])
                        

                        subplot(5,5,3)
                                                hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yHCNfield,'ShowText','on')
                        title('yHCN ')

                        pbaspect([1 1 1])
                        

                        subplot(5,5,4)
                                                hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yNH3field,'ShowText','on')
                        title('yNH3')

                        pbaspect([1 1 1])
                       

                        subplot(5,5,5)
                                                hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yH2field,'ShowText','on')
                        title('yH2')

                        pbaspect([1 1 1])
                       

                        disp('20% of optimization completed')
                        
                        
                         %% FEED Pressure PAIR
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            str(1).p=pressure(j);  %Calculate everything           
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
                            yH2field(i,j)=real(str(5).yH2);

                        end

                    end
                    
                         %normalize results
                     
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                     
                    subplot(5,5,6)
                     hold on 
                     scatter(persFCH4,persPressure,'red','x')
                    contour(FCH4,pressure,yieldfield,'ShowText','on')
                 

                    pbaspect([1 1 1])
                    
                    ylabel('Pressure')

                    subplot(5,5,7)
                       hold on 
                     scatter(persFCH4,persPressure,'red','x')
                    contour(FCH4,pressure,pricefield,'ShowText','on')
              

                    pbaspect([1 1 1])
               

                    subplot(5,5,8)
                       hold on 
                     scatter(persFCH4,persPressure,'red','x')
                    contour(FCH4,pressure,yHCNfield,'ShowText','on')
                

                    pbaspect([1 1 1])
           

                    subplot(5,5,9)
                       hold on 
                     scatter(persFCH4,persPressure,'red','x')
                    contour(FCH4,pressure,yNH3field,'ShowText','on')
            

                    pbaspect([1 1 1])
    

                    subplot(5,5,10)
                       hold on 
                     scatter(persFCH4,persPressure,'red','x')
                    contour(FCH4,pressure,yH2field,'ShowText','on')
         

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('Pressure')
                    
                    
                    disp('40% of optimization completed')
                        
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
                            yH2field(i,j)=real(str(5).yH2);

                        end

                    end
                    
                          %normalize results
                    
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                    
                    subplot(5,5,11)
                       hold on 
                     scatter(persFCH4,persnrow,'red','x')
                    contour(FCH4,nrow,yieldfield,'ShowText','on')
                 

                    pbaspect([1 1 1])
               
                    ylabel('PFR Segments')

                    subplot(5,5,12)
                      hold on 
                     scatter(persFCH4,persnrow,'red','x')
                    contour(FCH4,nrow,pricefield,'ShowText','on')
                 

                    pbaspect([1 1 1])
          

                    subplot(5,5,13)
                      hold on 
                     scatter(persFCH4,persnrow,'red','x')
                    contour(FCH4,nrow,yHCNfield,'ShowText','on')
                   

                    pbaspect([1 1 1])
  
                    subplot(5,5,14)
                      hold on 
                     scatter(persFCH4,persnrow,'red','x')
                    contour(FCH4,nrow,yNH3field,'ShowText','on')
                

                    pbaspect([1 1 1])


                    subplot(5,5,15)
                      hold on 
                     scatter(persFCH4,persnrow,'red','x')
                    contour(FCH4,nrow,yH2field,'ShowText','on')
         
                    pbaspect([1 1 1])
                    
                    disp('60% of optimization completed')
                    
                                        %% FEED T PAIR
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            str(1).T=temperature(j);   %Calculate everything           
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
                            yH2field(i,j)=real(str(5).yH2);

                        end

                    end
                    
                     %normalize results
                    
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
                    

                    subplot(5,5,16)
                    contour(FCH4,temperature,yieldfield,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('Temperature Reactor inlet')

                    subplot(5,5,17)
                    contour(FCH4,temperature,pricefield,'ShowText','on')
                   

                    pbaspect([1 1 1])
              
                

                    subplot(5,5,18)
                    contour(FCH4,temperature,yHCNfield,'ShowText','on')
                  

                    pbaspect([1 1 1])
               
                

                    subplot(5,5,19)
                    contour(FCH4,temperature,yNH3field,'ShowText','on')
                 

                    pbaspect([1 1 1])
                
     

                    subplot(5,5,20)
                    contour(FCH4,temperature,yH2field,'ShowText','on')
                   

                    pbaspect([1 1 1])
                   
                    disp('80% of optimization completed')
                    
                    
                    %% FEED HEATMED PAIR
                    switch unt(1).cocross
                        case 1
                    for i=1:n

                        for j=1:n
                            str(1).FCH4=FCH4(i);
                            str(4).G=hstream(j);  %Calculate everything           
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
                            yH2field(i,j)=real(str(5).yH2);

                        end

                    end
                    
                    
                    
                     %normalize results
                
                     pricefield=pricefield/(max(pricefield(:))-min(pricefield(:)));
              

                    subplot(5,5,21)
                    contour(FCH4,hstream,yieldfield,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    ylabel('Molar Flow Heating Medium')

                    subplot(5,5,22)
                    contour(FCH4,hstream,pricefield,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                

                    subplot(5,5,23)
                    contour(FCH4,hstream,yHCNfield,'ShowText','on')
                  

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                

                    subplot(5,5,24)
                    contour(FCH4,hstream,yNH3field,'ShowText','on')
                 

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
     

                    subplot(5,5,25)
                    contour(FCH4,hstream,yH2field,'ShowText','on')
                   

                    pbaspect([1 1 1])
                    xlabel('Feed CH4')
                    
                    disp('Optimization Completed')
                    
                        case 0
                            
                            disp('Cross current heating is active for modelling, flow of heating medium is constant and  was ignored for optimization')
                            disp('Optimization Completed')


                    end
end