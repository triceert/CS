
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

%%%%%%%%%


%%

function main(varargin) %give options for what to execute how

%% get all subfolders and containing functions
mfilename='main.m';
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


%% RUN MODE (PROVISORISCH)
ir=1;               %IDEAL REAL 0 ideal/1 real plant modelling
pr=1;               %pressure drop (0 without, 1 with)
cc=0;               %crosscurrent 0 cocurrent 1  heating of reactor   
sens=1;             %sensitivity analysis on or off (0/1) TIME CONSUMING (5min plus)



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
    %latex 
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot,'defaultTextInterpreter','latex');
    
%EVALUATE CALL MODE AND DO IT ACCORDINGLY
    switch nargin    
        case 0    %if no argument, just do everything   
            [cmp,unt,str]=loader(excelstring);
            [cmp,unt,str]=calculator(cmp,unt,str,ir,cc,sens,pr);
            evaluator(cmp,unt,str);
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
                [~,~,~]=calculator(cmp,unt,str,ir,cc,sens,pr);
            elseif  comp3==1    %argument eval, check if calc data here and eval
                cmpcalc = evalin('base', 'cmpcalc');
                untcalc = evalin('base', 'untcalc');
                strcalc = evalin('base', 'strcalc');
                evaluator(cmpcalc,untcalc,strcalc);
            elseif comp4==1        %argument calceval, check if data here and calc and eval
                cmp = evalin('base', 'cmp');
                unt = evalin('base', 'unt');
                str = evalin('base', 'str');
                [cmpcalc,untcalc,strcalc]=calculator(cmp,unt,str,ir,cc,sens,pr);
                evaluator(cmpcalc,untcalc,strcalc);
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
function [cmp, unt, str]=calculator(cmp,unt,str,ir,cc,sens,pr)
    time2=tic;
    cprintf('blue','Calculations started\n');
    

%% INITIAL GUESS VALUES FOR ONE CONSISTENT RUN/OPTIMIZATION

%COCROSS AND IDEAL REAL SWITCH

%MODELLING PARAMETERS
unt(1).cocross=cc;           % 0 cross 1 co 
unt(1).ideal_real=ir;        %ideal gas 0, peng robinson 1
unt(1).pre_drop=pr;
%THERMO
str(1).p=103325;   %feed pressure
str(1).T=700.08;       %feed temperature
str(2).T=1600;      %touter 
%FEEDS
str(1).FCH4=0.0123873; %absolute feed ch4 per single tube mol s-1   
str(1).ubsch=1.09935;  %überschuss NH3
str(4).G=0.01;            %Heating Medium Flow rate, only usefull when co-current flow
                                %in flow heating medium per fucking tube(ignored if cross heated)
 %Reactor length
unt(1).nrow=1;             %number reactr elements in row

%flow ratio
str(10).L=520; %liquid flow rate stream 10 

%% CALL Optimizor (Böser Strom und Zeitfresser!)
senspara=sens;                     %0 sensitivity analysis off   1 on
[cmp,unt,str]=optimizor(cmp,unt,str,senspara);

 %% MAKE ONE CONSISTENT RUN AND PLOT   

        [cmp,unt,str]=reactorcalculator(cmp,unt,str,1);       %plotparamter 1
        [cmp,unt,str]=nh3absorber(cmp,unt,str);       
        [cmp,unt,str] = hcnabsorber(cmp,unt,str);         
        [cmp,unt,str]=hcn_distillation(cmp,unt,str);
        [cmp,unt,str] = hcnabsorber(cmp,unt,str); 
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

function evaluator(cmpin, untin, strin)
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
fprintf('Column diameter [m]: D = %g\n', untin(2).dia);
fprintf('Flow rate ratios: L/G = %g\n', untin(2).ratio);
fprintf('Volumetric fraction of h2so4: phi = %g\n', untin(2).frac_h2so4);
fprintf('NH3 Column CAPEX : Capex = %g\n', untin(2).capex);
fprintf('NH3 Column OPEX : Opex = %g\n', untin(2).opex);
fprintf('NH3 Column Temperature: T = %g\n', strin(7).T);
        

cprintf('Blue','HCN Adsorber:\n')       
fprintf('Number of theoretical units: NTU = %g\n', untin(3).ntu);
fprintf('Height of theoretical units: HTU = %g\n',  untin(3).htu);
fprintf('Column Height [m]: H = %g\n', untin(3).h);
fprintf('Column diameter [m]: D = %g\n', untin(3).dia);
fprintf('Flow rate ratios: L/G = %g\n', untin(3).ratio);
fprintf('HCN Absorber CAPEX : Capex = %g\n',  untin(3).capex);
fprintf('HCN Absorber OPEX : Opex = %g\n', untin(3).opex);
fprintf('HCN Absorber Outlet tempeature [K] : Temp = %g\n', strin(9).T);
fprintf('HCN Absorber Outlet fraction of HCN : x_HCN = %g\n', strin(9).xHCN);
fprintf('Inlet molar flow of liquid stream : L_in = %g\n', strin(10).L);
        
cprintf('Blue','Distillation:\n')
fprintf('Number of stages []: H = %g\n', untin(4).N_S);
fprintf('Column Height [m]: H = %g\n', untin(4).h);
fprintf('Column Diameter [m]: d = %g\n', 2*untin(4).rad);
fprintf('Condenser heat duty [W]: Q_C = %g\n', untin(4).Q_cond);
fprintf('Reboiler heat duty [W]: Q_C = %g\n', untin(4).Q_reb);
fprintf('HX bef. dist. heat duty [W]: Q_C = %g\n', untin(4).Q_HX_before_dist);
fprintf('Condenser area [m2]: A_C = %g\n', untin(4).A_cond);
fprintf('Reboiler area [m2]: A_R = %g\n', untin(4).A_reb);
fprintf('HX bef. dist. area [m2]: A_HX_bef_dist = %g\n', untin(4).A_HX_before_dist);
fprintf('optimized rate [mol/s]: optimzed flow rate = %g\n', strin(10).L);
fprintf('Dist Capex/$ = %g\n', untin(4).capex);
fprintf('Dist Opex/$/y = %g\n', untin(4).opex);

 %assign for outputs UNITS
 unttable=untin;
 fields={'As','rad','Xopt','Aq','deltaHHCN','Reout','ideal_real','deltaHNH','cocross','idgc','GHV','En','totex','price_processwater','price_coolingwater','price_natural_gas','nrow','N_tubes',...
     'N_tubes_aside','HCNneeded','a','Q_in','htu','rho_mix_in','Q_tot','conv','toex','tspan','lspan','epsilon','price_tube','price_vapor','deltaWall'};
 unttable = rmfield(unttable,fields);
 unttable(5).cocross=untin(1).cocross;
 unttable(5).idealreal=untin(1).ideal_real;
 unttable(1).h=untin(1).nrow*2; %reassign new reactor length
 unttable(7)=[];
 
UNT_Table=struct2table(unttable);


%assign Streams
strtable=strin;
fields={'FCH4','ubsch','xCH4','xNH3','xH2','T','yAS','yH2SO4','xN2'};
strtable = rmfield(strtable,fields);
strtable(13)=[];
STR_Table=struct2table(strtable);
 
writetable(UNT_Table,'Run.xls','Sheet',1)
writetable(STR_Table,'Run.xls','Sheet',2)
        

    cprintf('_comment','Plots and table generation successfull \n');
    cprintf('_comment','Plots saved as PDFs \n');
    cprintf('_comment','Run Excel sheet generated \n');
end









