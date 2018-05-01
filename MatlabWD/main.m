
%% MAIN FUNCTION CALLER
%Function Description: A function that calls and calculates everything
%Arguments:
%      Takes various number of arguments;
%       NO ARGUMENT: default execution, loads calculates and evaluates
%       everything
%      ONE ARGUMENT:
%           'load'      just load data
%           'calc'      just calculate, only works if data loaded
%           'eval'      just evaluate, only works if calculations made
%           'calceval'  calc and evaluate, only works if data loaded
%Outputs:
%      NIL

function main(varargin) %give options for what to execute how
%INIT
    %vars
    excelstring='Data.xlsx'; %Name of underlng master excel file
    %set begin timer
    time1=tic;    
    %clean up
    clear vars
    close all
    clc
    cprintf('_Blue','Execution started.\n');%say what we did here 
    
%EVALUATE CALL MODE AND DO IT ACCORDINGLY
    switch nargin    
        case 0    %if no argument, just do everything   
            [cmp,unt,str]=loader(excelstring);
            [cmp,unt,str]=calculator(cmp,unt,str);
            [~,~]=evaluator(cmp,unt,str);
        case 1    %if one argument, check what argument
            comp1=strcmp(varargin{1},'load');
            comp2=strcmp(varargin{1},'calc');
            comp3=strcmp(varargin{1},'eval');
            comp4=strcmp(varargin{1},'calceval');
            if comp1==1 %argument load, just load
                [~,~,~]=loader(excelstring);
            elseif  comp2==1 %argument calc, check if data available and calc
                cmp = evalin('base', 'cmp');
                unt = evalin('base', 'unt');
                str = evalin('base', 'str');
                [~,~,~]=calculator(cmp,unt,str);
            elseif  comp3==1    %argument eval, check if calc data here and eval
                cmpcalc = evalin('base', 'cmpcalc');
                untcalc = evalin('base', 'untcalc');
                strcalc = evalin('base', 'strcalc');
                [~,~]=evaluator(cmpcalc,untcalc,strcalc);
            elseif comp4==1        %argument calceval, check if data here and calc and eval
                cmp = evalin('base', 'cmp');
                unt = evalin('base', 'unt');
                str = evalin('base', 'str');
                [cmpcalc,untcalc,strcalc]=calculator(cmp,unt,str);
                [~,~,~]=evaluator(cmpcalc,untcalc,strcalc);
            else
                error(...
                    'No valid function argument in mainexec. For default call without argument.')
            end         
        otherwise
            error('No valid number of function arguments in mainexec. For default call without argument.')
    end        

%TERMINATE
    elapse=toc(time1);
    cprintf('_comment','Main execution terminated normally in a total of\f');
    cprintf('_comment','%g',elapse);
    cprintf('_comment','\f seconds.\n');%say what we did here&timer
end  

%% %%%%%%%%%%%%%%%
%NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%

%% LOADER
%Function Description: A function that calls and calculates everything
%Arguments:
%      string    string to the excel files    
%Outputs:
%      cmp unt   structs containing importet compound and unit operations
%      data
function [cmp,unt,str]=loader(string)
   [cmp,unt,str]=dataopener(string); %load data from file
   
   assignin('base','cmp',cmp)    %assign for everyone
   assignin('base','unt',unt) 
   assignin('base','str',str) 
end

%% Calculator
%Function Description: A function that calls and calculates everything
%Arguments:
%      cmpin,untin      compounds and unit operations data from anywhere
%Outputs:
%      cmpout,untout    compunds and unit operations data to anywhere after
%      calculation
function [cmp, unt, str]=calculator(cmp,unt,str)
    time2=tic;
    cprintf('blue','Calculations started\n');
    

%% INITIAL GUESS VALUES FOR ONE CONSISTENT RUN/OPTIMIZATION

%COCROSS AND IDEAL REAL SWITCH

%MODELLING PARAMETERS
unt(1).cocross=1;           % 0 cross 1 co 
unt(1).ideal_real=1;        %ideal gas 0, peng robinson 1
%THERMO
str(1).p=103325;   %feed pressure
str(1).T=699.999;       %feed temperature
str(2).T=1600;      %touter 
%FEEDS
str(1).FCH4=0.0776669; %absolute feed ch4 per single tube mol s-1   
str(1).ubsch=1.05;  %überschuss NH3
str(4).G=0.09;            %Heating Medium Flow rate, only usefull when co-current flow
                                %in flow heating medium per fucking tube(ignored if cross heated)
 %Reactor length
 unt(1).nrow=1;             %number reactr elements in row

%% CALL Optimizor (Böser Strom und Zeitfresser!)
senspara=1;                     %0 sensitivity analysis off   1 on
[cmp,unt,str]=optimizor(cmp,unt,str,senspara);

 %% MAKE ONE CONSISTENT RUN AND PLOT   

        [cmp,unt,str]=reactorcalculator(cmp,unt,str,1);       %plotparamter 1
        [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
        [cmp,unt,str] = hcnideal(cmp,unt,str);       
        [cmp,unt,str]=hcn_distillation(cmp,unt,str);       
        [cmp,unt,str]=pricecalculator(cmp,unt,str);
        
        
      %assign for everyone
        cmpout=cmp;
        strout=str;
        untout=unt;
        
   
        assignin('base','cmpcalc',cmpout)   
        assignin('base','untcalc',untout) 
        assignin('base','strcalc',strout) 

    %say what we did here
    elapse=toc(time2);
    cprintf('_comment',...
        'Stuff was calculated successfully. ...Pure calculation time was\n');
    cprintf('_comment','%g',elapse);
    cprintf('_comment','\f seconds.\n')  
end

%% EVALUATOR
%Function Description: A general function for plots and output generations
%Arguments:
%      mpin,untin      compounds and unit operations data from anywhere
%Outputs:
%      plots, table    plots and tables for export or whatever

function [tables, plots]=evaluator(cmpin, untin, strin)
    cprintf('blue','Begin to plot and generate export files\n');
    
 cprintf('Blue','Main run model output:\n')     
    
    

    

    
 

    

cprintf('Blue','Plant:\n')    

    %MODEL RUN unt(1).co unt(1).ideal_real
    fprintf('break even $/kg =%g',untin(5).price)
    fprintf('$\n')
    
cprintf('Blue','Reactor:\n')     

fprintf('Feed STream CH4 [mol/s] = %g\n', (strin(1).G*strin(1).yCH4));
fprintf('Feed STream CH4 per tube: [mol/s] = %g\n', (strin(1).FCH4));
fprintf('OutStream HCN [mol/s] = %g\n', (strin(5).G*strin(5).yHCN));
fprintf('InStream Heatmed [mol/s] = %g\n', (strin(4).G));
fprintf('yield HCN resp, CH4: FCH4 [] = %g\n', untin(1).yield);
fprintf('NH3 Excess ratio [] = %g\n', strin(1).ubsch);
fprintf('Conversion CH4 [] = %g\n', untin(1).conv);
fprintf('yCH4 outlet [] = %g\n', strin(5).yCH4);
fprintf('yNH3 outlet [] = %g\n', strin(5).yNH3);
fprintf('yHCN outlet [] = %g\n', strin(5).yHCN);
fprintf('yH2 outlet [] = %g\n', strin(5).yH2);
fprintf('yN2 outlet [] = %g\n', strin(5).yN2);
fprintf('Temperature inlet [] = %g\n', strin(1).T);
fprintf('pressure inlet= %g\n', strin(1).p);


fprintf('Reactor Length [m] = %g\n', untin(1).nrow*2);
fprintf('Tubes needed = %g\n', untin(1).N_tubes);
fprintf('Capex = %g\n', untin(1).capex);
fprintf('Oapex = %g\n', untin(1).opex);


    
        
     
cprintf('Blue','NH3 Adsorber:\n') 
fprintf('Number of theoretical units: NTU = %g\n', untin(2).ntu);
fprintf('Height of theoretical units: HTU = %g\n', untin(2).htu);
fprintf('Column Height [m]: H = %g\n', untin(2).h);

fprintf('NH3 Column CAPEX [Mio. US$]: Capex = %g\n', untin(2).capex);
fprintf('NH3 Column OPEX [Mio. US$]: Opex = %g\n', untin(2).opex);
        

cprintf('Blue','HCN Adsorber:\n')       
fprintf('Number of theoretical units: NTU = %g\n', untin(3).ntu);
fprintf('Height of theoretical units: HTU = %g\n',  untin(3).htu);
fprintf('Column Height [m]: H = %g\n', untin(3).h);
fprintf('HCN Column CAPEX [Mio. US$]: Capex = %g\n',  untin(3).capex);
fprintf('HCN Column OPEX [Mio. US$]: Opex = %g\n', untin(3).opex);

        
cprintf('Blue','Distillation:\n')          
 fprintf('Column Height [m]: H = %g\n', untin(4).h);


        
        %Call Plotter Function
        harry_plotter
        
        
        
        %tables= call tex table maker
        %call some economic evaluator
        tables=NaN; %dummy
        plots=NaN;  %dummy
        assignin('base','plots',plots)   
        assignin('base','textable',tables) 

    cprintf('_comment','Plots and table generation successfull\n');
end









