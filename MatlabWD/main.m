
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
function [cmpout, untout, strout]=calculator(cmpin,untin,strin)
    time2=tic;
    cprintf('blue','Calculations started\n');
    
    %calculate different stuff
    
            
%idealreal 0 ideal 1 real(Peng robinson)


% idealreal=1; %best for real ohne t profile
% Pressure=101325;
% FeedCH4=   0.0010;   %varieren zwischen 0.001 und 0.1
% uberschuss=1.05;
% Tfeed=700;
% Touter=1600;
% pfrseries= 3;       %varieren zwischen 1 und 10

%set parameters BEST FOR REAL ohne T profile


% idealreal=0; %best for ideal
% Pressure=101325;
% FeedCH4=  0.01;
% uberschuss=1.05;
% Tfeed=700;
% Touter=1600;
% pfrseries=  1;


%% VALUES FOR ONE CONSISTENT RUN


%COCROSS AND IDEAL REAL SWITCH





idealreal=1;
cocross=0;           % 0 cross 1 co
strin(1).ideal_real=idealreal;
strin(1).cocross=cocross;

strin(4).G=0.0010;      %flow heating medium per fucking tube(ignored if cross heated)

strin(1).p=101325;   %feed pressure
strin(1).T=700;       %feed temperature
strin(2).T=600;      %touter 

strin(1).FCH4=0.0010; %absolute feed ch4 per single tube mol s-1   
strin(1).ubsch=1.05;  %überschuss NH3

untin(1).nrow=1;      %number reactr elements in row
 
[cmp,unt,str]=reactorcalculator(cmpin,untin,strin,1);%test
%CALL Optimizor
%optimizor(cmpin,untin,strin,idealreal,cocross);
                   




 %% MAKE ONE CONSISTENT RUN AND PLOT   
    
    
    
        
        %[cmp,unt,str]=reactorcalculator(cmpin,untin,strin,1);       %plotparamter 1
        [cmp,unt,str]=NH3_absorber_ideal(cmp,unt,str);       
        [cmp,unt,str] = hcnabsorption2(cmp,unt,str);       
        [cmp,unt,str]=hcn_distillation(cmp,unt,str); 
        
        [unt,str]=OPEX_reactor(cmp,unt,str);
        [unt]=CAPEX_reactor(unt);
        [unt]=TOTEX_reactor(unt);
        break_even_price=pricecalculator(unt,cmp)
        

        
        cmpout=cmp;
        strout=str;
        untout=unt;
        
     %assign for everyone
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

function [tables, plots]=evaluator(cmpin, untin,strin)
    cprintf('blue','Begin to plot and generate export files\n');
    
    
%         untout(2).htu = HTU;    %height theoretical units NH3 ADSORBER
%         untout(2).ntu = NTU;    %number theoritac units
%         untout(2).h = Z;        %height NH3 Absorber
%         untout(2).V = V_column; %VOlume
%         untout(2).capex = capex;
%         untout(2).opex = opex_tot2;
%         
%         untout(3).h = h;       %height             %HCN ABSORBTION
%         untout(3).capex = CAPEX_mil; %capex
%         untout(3).opex = opex_tot;   %opex
%         untout(3).ntu = NTU;            %theoretical number units   
%         untout(3).htu = HTU;            %height tu
%         untout(3).V = V_column;         %volume

        
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









