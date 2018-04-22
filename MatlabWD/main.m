
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
            [cmp,unt]=loader(excelstring);
            [cmp,unt]=calculator(cmp,unt);
            [~,~]=evaluator(cmp,unt);
        case 1    %if one argument, check what argument
            comp1=strcmp(varargin{1},'load');
            comp2=strcmp(varargin{1},'calc');
            comp3=strcmp(varargin{1},'eval');
            comp4=strcmp(varargin{1},'calceval');
            if comp1==1 %argument load, just load
                [~,~]=loader(excelstring);
                disp('load')
            elseif  comp2==1 %argument calc, check if data available and calc
                cmp = evalin('base', 'cmp');
                unt = evalin('base', 'unt');
                [~,~]=calculator(cmp,unt);
                disp('calc')
            elseif  comp3==1    %argument eval, check if calc data here and eval
                cmpcalc = evalin('base', 'cmpcalc');
                untcalc = evalin('base', 'untcalc');
                [~,~]=evaluator(cmpcalc,untcalc);
                disp('eval')
            elseif comp4==1        %argument calceval, check if data here and calc and eval
                cmp = evalin('base', 'cmp');
                unt = evalin('base', 'unt');
                [cmpcalc,untcalc]=calculator(cmp,unt);
                [~,~]=evaluator(cmpcalc,untcalc);
                disp('calceval')
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
function [cmp,unt]=loader(string)
   [cmp,unt]=dataopener(string); %load data from file
   
   assignin('base','cmp',cmp)    %assign for everyone
   assignin('base','unt',unt) 
end

%% Calculator
%Function Description: A function that calls and calculates everything
%Arguments:
%      cmpin,untin      compounds and unit operations data from anywhere
%Outputs:
%      cmpout,untout    compunds and unit operations data to anywhere after
%      calculation
function [cmpout, untout]=calculator(cmpin,untin)
    time2=tic;
    cprintf('blue','Calculations started\n');
    
    %calculate different stuff
        %Calc Function 1
        %Calc Function 2
        %Calc Function 3
        cmpout=cmpin; %dummy
        untout=untin; %dummy
        
     %assign for everyone
        assignin('base','cmpcalc',cmpout)   
        assignin('base','untcalc',untout) 

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

function [tables, plots]=evaluator(cmpin, untin)
    cprintf('blue','Begin to plot and generate export files\n');
    
        %plots=Call Plotter Function
        %tables= call tex table maker
        %call some economic evaluator
        tables=NaN; %dummy
        plots=NaN;  %dummy
        assignin('base','plots',plots)   
        assignin('base','textable',tables) 

    cprintf('_comment','Plots and table generation successfull\n');
end





%%% commented out, eventually usefull
%      if evalin('base','exist(''cmp'',''var'')')==1        
%      else
%          warning('Data for calculations do not exist')
%      end
%     




