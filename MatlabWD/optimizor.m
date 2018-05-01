function [cmp,unt,str]=optimizor(cmp,unt,str)


 %% OPITMIZORRR

    %initialize
    disp('Optimization started, dependent on the performance of your \r\n computer this can take up to 5 minutes.\r\n It is a good time to grab a coffe now.')

    n=4; %optim grid size n*n..dont push it to far up
    
    %get persistent/guess variables
    persFCH4=str(1).FCH4;
    persRatio=str(1).ubsch;
    persPressure=str(1).p;
    persTemp=str(1).T;
    persnrow=unt(1).nrow;   
    persFlow=str(4).G;


    
%% OPTIMIZE    
    
    function price=optimfunc(x,cmp,unt,str)  %Price Optimization Objective Function
        %get initial guess
          str(1).FCH4=x(1);
          str(1).ubsch=x(2);
          str(1).p=x(3);
          str(1).T=x(4);
          unt(1).nrow=x(5);       
          str(4).G=x(6);
          %calculate everything
                                [cmp,unt,str]=reactorcalculator(cmp,unt,str,0);   %plotparameter 0
                                [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
                                [cmp,unt,str] = hcnideal(cmp,unt,str);       
                                [cmp,unt,str]=hcn_distillation(cmp,unt,str); 
                                [unt,str]=OPEX_reactor(cmp,unt,str);
                                [unt]=CAPEX_reactor(unt);
                                [unt]=TOTEX_reactor(unt);
                                 price=pricecalculator(unt,cmp);                             
    end


   %Find Optimum
     %fmincon initialize
     optimhandle=@(x)optimfunc(x,cmp,unt,str);
     options=optimset('Display','off');
      x0 = [persFCH4,persRatio,persPressure,persTemp,persnrow,persFlow];
      lb = [1e-4,0.9,85000,650,1,1e-4];
      ub = [1e-1,1.1,500000,750,15,1];
      disp('Fmincon startup complete.')
      x = fmincon(optimhandle,x0,[],[],[],[],lb,ub,[],options); %solve
      disp('20% of optimization completed')
            persFCH4=x(1);
            persRatio=x(2);
            persPressure=x(3);
            persTemp=x(4);
            persnrow=x(5);
            persnrow=round(persnrow); %can only be discrete
            persFlow=x(6);
            
                str(1).FCH4=persFCH4;
                str(1).ubsch=persRatio;
                str(1).p=persPressure;
                str(1).T=persTemp;
                unt(1).nrow=persnrow;
                str(4).G=persFlow;
      disp('Fmincon successfully found optimal values.\n Continuing with sensitivity analysis.')
  
                
 %% SENSITIVITY ANALYSIS STARTUP  
    %a=0.2; %low up ratio %
    %lb=x-x*a; %lower bound for sensitivity
    %ub=x+x*a;
    strprov=str;
    untprov=unt;
    cmpprov=cmp; 
    
    
    %Assign range for sensitivity analysis
    FCH4 = linspace(lb(1),ub(1),n);    %FEEDRANGE
    ubsch=linspace(lb(2),ub(2),n);       %EXCESS NH3 RANGE
    pressure=linspace(lb(3),ub(3),n);         %Pressure
    temperature=linspace(lb(4),ub(5),n);     %Temoerature
        %if persnrow<4
            %lb(5)=1;
        %else
            %lb(5)=persnrow-3;
        %end  
        %ub(5)=persnrow+3;
    nrow=linspace(lb(5),ub(5),n);                %NUMBER OF PFRs in row range
    hstrproveam=linspace(lb(6),ub(6),n);           %HEAT strprovEAM MOLAR FLOW
  

  
     

    %preallocate memory for different sensitivity parameters
    yieldfield=zeros(numel(n,n));
    conversionfield=zeros(numel(n,n));
    yNH3field=zeros(numel(n,n));
    yH2field=zeros(numel(n,n));
    pricefield=zeros(numel(n,n));
    yHCNfield=zeros(numel(n,n));
    
      

            %% FEED EXCESS PAIR
                        for i=1:n
                            for j=1:n
                                strprov(1).FCH4=FCH4(i);
                                strprov(1).ubsch=ubsch(j);  %Calculate everything           
                                        [cmpprov,untprov,strprov]=reactorcalculator(cmpprov,untprov,strprov,0);   %plotparameter 0
                                        [cmpprov,untprov,strprov]=NH3_absorber_ideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov] = hcnideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov]=hcn_distillation(cmpprov,untprov,strprov); 
                                        [untprov,strprov]=OPEX_reactor(cmpprov,untprov,strprov);
                                        [untprov]=CAPEX_reactor(untprov);
                                        [untprov]=TOTEX_reactor(untprov);
                                pricefield(i,j)=pricecalculator(untprov,cmpprov);       %assign into mesh
                                yieldfield(i,j)=real(untprov(1).yield);                                                                                                                    
                                yNH3field(i,j)=real(strprov(5).yNH3);
                                yH2field(i,j)=real(strprov(5).yH2);
                                yHCNfield(i,j)=real(strprov(5).yHCN);
                            end
                        end
                     
                    % normalize results   
                     pricefield=(pricefield-min(pricefield(:)))/(max(pricefield(:))-min(pricefield(:)));
                      yieldfield=(yieldfield-min(yieldfield(:)))/(max(yieldfield(:))-min(yieldfield(:)));
                      yNH3field=(yNH3field-min(yNH3field(:)))/(max(yNH3field(:))-min(yNH3field(:)));
                      yHCNfield=(yHCNfield-min(yHCNfield(:)))/(max(yHCNfield(:))-min(yHCNfield(:)));
                      yH2field=(yH2field-min(yH2field(:)))/(max(yH2field(:))-min(yH2field(:)));
                     
                    
%plots
                    figure(1)
                        subplot(4,2,1)
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yieldfield,'ShowText','on')
                        title('Yield  resp. to CH4')
                        pbaspect([1 1 1])                    
                        ylabel('Excess NH3')
                        
                        subplot(4,2,2)                     
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,pricefield,'ShowText','on')
                        title('Break even Price $/kg')
                        pbaspect([1 1 1])                        

       
                        
                        figure(2)
                        subplot(4,3,1)
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yNH3field,'ShowText','on')
                        title('yNH3')
                        pbaspect([1 1 1])  
                        ylabel('Excess NH3')

                        subplot(4,3,2)
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yH2field,'ShowText','on')
                        title('yH2')
                        pbaspect([1 1 1])
                        
                        subplot(4,3,3)
                        hold on 
                        scatter(persFCH4,persRatio,'red','x')
                        contour(FCH4,ubsch,yHCNfield,'ShowText','on')
                        title('yHCN')
                        pbaspect([1 1 1])
                       

                        disp('30% of optimization completed')
   %% Feed Pressure Pair                     
                        for i=1:n
                            for j=1:n
                                strprov(1).FCH4=FCH4(i);
                                strprov(1).pressure=ubsch(j);  %Calculate everything           
                                        [cmpprov,untprov,strprov]=reactorcalculator(cmpprov,untprov,strprov,0);   %plotparameter 0
                                        [cmpprov,untprov,strprov]=NH3_absorber_ideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov] = hcnideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov]=hcn_distillation(cmpprov,untprov,strprov); 
                                        [untprov,strprov]=OPEX_reactor(cmpprov,untprov,strprov);
                                        [untprov]=CAPEX_reactor(untprov);
                                        [untprov]=TOTEX_reactor(untprov);
                                pricefield(i,j)=pricecalculator(untprov,cmpprov);       %assign into mesh
                                yieldfield(i,j)=real(untprov(1).yield);                                                                                     
                                conversionfield(i,j)=real(untprov(1).conv);
                                yNH3field(i,j)=real(strprov(5).yNH3);
                                yH2field(i,j)=real(strprov(5).yH2);
                                yHCNfield(i,j)=real(strprov(5).yHCN);
                            end
                        end
                     
                    % normalize results   
                     pricefield=(pricefield-min(pricefield(:)))/(max(pricefield(:))-min(pricefield(:)));
                      yieldfield=(yieldfield-min(yieldfield(:)))/(max(yieldfield(:))-min(yieldfield(:)));
                      yNH3field=(yNH3field-min(yNH3field(:)))/(max(yNH3field(:))-min(yNH3field(:)));
                      yHCNfield=(yHCNfield-min(yHCNfield(:)))/(max(yHCNfield(:))-min(yHCNfield(:)));
                      yH2field=(yH2field-min(yH2field(:)))/(max(yH2field(:))-min(yH2field(:)));
                     
                    
%plots
                    figure(1)
                        subplot(4,2,3)
                        hold on 
                        scatter(persFCH4,persPressure,'red','x')
                        contour(FCH4,pressure,yieldfield,'ShowText','on')
                        
                        pbaspect([1 1 1])                    
                        ylabel('Pressure')
                        
                        subplot(4,2,4)                     
                        hold on 
                        scatter(persFCH4,persPressure,'red','x')
                        contour(FCH4,pressure,pricefield,'ShowText','on')
                        
                        pbaspect([1 1 1])                        

       
                        
                        figure(2)
                        subplot(4,3,4)
                        hold on 
                        scatter(persFCH4,persPressure,'red','x')
                        contour(FCH4,pressure,yNH3field,'ShowText','on')
                        ylabel('Pressure')
                        pbaspect([1 1 1])                     

                        subplot(4,3,5)
                        hold on 
                        scatter(persFCH4,persPressure,'red','x')
                        contour(FCH4,pressure,yH2field,'ShowText','on')
              
                        pbaspect([1 1 1])
                        
                        subplot(4,3,6)
                        hold on 
                        scatter(persFCH4,persPressure,'red','x')
                        contour(FCH4,pressure,yHCNfield,'ShowText','on')
            
                        pbaspect([1 1 1])
                       

                        
              
                    
                    %% FEED HEATMED PAIR
                    switch untprov(1).cocross
                        case 1
                            
                                           
                        for i=1:n
                            for j=1:n
                                strprov(1).FCH4=FCH4(i);
                                strprov(4).G=hstrproveam(j);  %Calculate everything           
                                        [cmpprov,untprov,strprov]=reactorcalculator(cmpprov,untprov,strprov,0);   %plotparameter 0
                                        [cmpprov,untprov,strprov]=NH3_absorber_ideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov] = hcnideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov]=hcn_distillation(cmpprov,untprov,strprov); 
                                        [untprov,strprov]=OPEX_reactor(cmpprov,untprov,strprov);
                                        [untprov]=CAPEX_reactor(untprov);
                                        [untprov]=TOTEX_reactor(untprov);
                                pricefield(i,j)=pricecalculator(untprov,cmpprov);       %assign into mesh
                                yieldfield(i,j)=real(untprov(1).yield);                                                                                     
                                conversionfield(i,j)=real(untprov(1).conv);
                                yNH3field(i,j)=real(strprov(5).yNH3);
                                yH2field(i,j)=real(strprov(5).yH2);
                                yHCNfield(i,j)=real(strprov(5).yHCN);
                            end
                        end
                     
                    % normalize results   
                     pricefield=(pricefield-min(pricefield(:)))/(max(pricefield(:))-min(pricefield(:)));
                      yieldfield=(yieldfield-min(yieldfield(:)))/(max(yieldfield(:))-min(yieldfield(:)));
                      yNH3field=(yNH3field-min(yNH3field(:)))/(max(yNH3field(:))-min(yNH3field(:)));
                      yHCNfield=(yHCNfield-min(yHCNfield(:)))/(max(yHCNfield(:))-min(yHCNfield(:)));
                      yH2field=(yH2field-min(yH2field(:)))/(max(yH2field(:))-min(yH2field(:)));
                     
                     
                    
%plots
                    figure(1)
                        subplot(4,2,7)
                        hold on 
                        scatter(persFCH4,persFlow,'red','x')
                        contour(FCH4,hstrproveam,yieldfield,'ShowText','on')
                    
                        pbaspect([1 1 1])                    
                        ylabel('Flow heating medium p. tube.')
                        
                        subplot(4,2,8)                     
                        hold on 
                        scatter(persFCH4,persFlow,'red','x')
                        contour(FCH4,hstrproveam,pricefield,'ShowText','on')
          
                        pbaspect([1 1 1])                        


                        
                        figure(2)
                        subplot(4,3,7)
                        hold on 
                        ylabel('Flow heating medium p. tube.')
                        scatter(persFCH4,persFlow,'red','x')
                        contour(FCH4,hstrproveam,yNH3field,'ShowText','on')
           
                        pbaspect([1 1 1])                     

                        subplot(4,3,8)
                        hold on 
                        scatter(persFCH4,persFlow,'red','x')
                        contour(FCH4,hstrproveam,yH2field,'ShowText','on')
   
                        pbaspect([1 1 1])
                        
                        subplot(4,3,9)
                        hold on 
                        scatter(persFCH4,persFlow,'red','x')
                        contour(FCH4,hstrproveam,yHCNfield,'ShowText','on')
          
                        pbaspect([1 1 1])
                    
                        disp('80% of optimization completed')
                        case 0
                            
                            disp('Cross current heating is active for modelling, flow of heating medium is constant and  was ignored for optimization')
                        disp('80% of optimization completed')   


                    end
                    
                    
%% nrow feed pair


          
                            
                                           
                        for i=1:n
                            for j=1:n
                                strprov(1).FCH4=FCH4(i);
                                untprov(1).nrow =nrow(j);  %Calculate everything           
                                        [cmpprov,untprov,strprov]=reactorcalculator(cmpprov,untprov,strprov,0);   %plotparameter 0
                                        [cmpprov,untprov,strprov]=NH3_absorber_ideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov] = hcnideal(cmpprov,untprov,strprov);       
                                        [cmpprov,untprov,strprov]=hcn_distillation(cmpprov,untprov,strprov); 
                                        [untprov,strprov]=OPEX_reactor(cmpprov,untprov,strprov);
                                        [untprov]=CAPEX_reactor(untprov);
                                        [untprov]=TOTEX_reactor(untprov);
                                pricefield(i,j)=pricecalculator(untprov,cmpprov);       %assign into mesh
                                yieldfield(i,j)=real(untprov(1).yield);                                                                                     
                                conversionfield(i,j)=real(untprov(1).conv);
                                yNH3field(i,j)=real(strprov(5).yNH3);
                                yH2field(i,j)=real(strprov(5).yH2);
                                yHCNfield(i,j)=real(strprov(5).yHCN);
                            end
                        end
                    % normalize results   
                     pricefield=(pricefield-min(pricefield(:)))/(max(pricefield(:))-min(pricefield(:)));
                      yieldfield=(yieldfield-min(yieldfield(:)))/(max(yieldfield(:))-min(yieldfield(:)));
                      yNH3field=(yNH3field-min(yNH3field(:)))/(max(yNH3field(:))-min(yNH3field(:)));
                      yHCNfield=(yHCNfield-min(yHCNfield(:)))/(max(yHCNfield(:))-min(yHCNfield(:)));
                      yH2field=(yH2field-min(yH2field(:)))/(max(yH2field(:))-min(yH2field(:)));
                     
                 switch untprov(1).cocross
                        case 0
                     
%plots
                    figure(1)
                        subplot(4,2,5)
                        hold on 
                        scatter(persFCH4,persnrow,'red','x')
                        contour(FCH4,nrow,yieldfield,'ShowText','on')
                        xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])                    
                        ylabel('n segements in row')
                        
                        subplot(4,2,6)                     
                        hold on 
                        scatter(persFCH4,persnrow ,'red','x')
                        contour(FCH4,nrow,pricefield,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])                        

                        
                        figure(2)
                        subplot(4,3,7)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yNH3field,'ShowText','on')
                         xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])    
                          ylabel('n segements in row')

                        subplot(4,3,8)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yH2field,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])
                        
                        subplot(4,3,9)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yHCNfield,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])
                     
                        
             case 1
                                                figure(1)
                        subplot(4,2,7)
                        hold on 
                        scatter(persFCH4,persnrow,'red','x')
                        contour(FCH4,nrow,yieldfield,'ShowText','on')
                        xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])                    
                        ylabel('n segements in row')
                        
                        subplot(4,2,8)                     
                        hold on 
                        scatter(persFCH4,persnrow ,'red','x')
                        contour(FCH4,nrow,pricefield,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])                        


                        
                        figure(2)
                        subplot(4,3,10)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yNH3field,'ShowText','on')
                         xlabel('Feed CH4 p.t.')
                         ylabel('n segements in row')
                        pbaspect([1 1 1])                     

                        subplot(4,3,11)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yH2field,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])
                        
                        subplot(4,3,12)
                        hold on 
                        scatter(persFCH4, persnrow,'red','x')
                        contour(FCH4,nrow,yHCNfield,'ShowText','on')
                          xlabel('Feed CH4 p.t.')
                        pbaspect([1 1 1])
                         
                                
                                
                 end
                 
                  disp('Optimization Completed')
end