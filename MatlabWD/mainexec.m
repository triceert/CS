
%% MAIN FUNCTION CALLER
%Function Description: A function that calls and calculates everything
%Arguments:
%      NIL      Does nothing so far
%Outputs:
%      a        control variable
%%

function mainexec %give options for what to execute how


%INIT
    %set begin timer
    time1=tic;
    
    %clean up
    clear vars
    close all
    clc
   
    
    cprintf('_blue','Workspace is tidy now. Execution started.\n');%say what we did here


    %call data loader and write inital compound and unit
    %operation info to local workspace of mainexec
    [cmp,uts]=dataopener('Data.xlsx');
    
    assignin('caller','compounds',cmp)%write structs to workspace
    assignin('caller','unit_operations',uts)%write structs to workspace
   
    cprintf('_comment','Data loaded successfully\n');%say what we did here

%CALCULATE STUFF
    time2=tic;
    cprintf('blue','Calculations started\n');
    
    %calculate different stuff


    %say what we did here
    elapse=toc(time2);
    cprintf('_comment','Stuff was calculated successfully. Pure calculation time was\n');
    cprintf('_comment','%g',elapse);
    cprintf('_comment','\f seconds.\n')


%EVALUATE STUFF
    cprintf('blue','Begin to plot and generate files\n');
    %call data closer&document writer




    cprintf('_comment','Plots and table generation successfull\n');


%TERMINATE

    elapse=toc(time1);
    cprintf('_comment','Main execution terminated normally in a total of\f');
    cprintf('_comment','%g',elapse);
    cprintf('_comment','\f seconds.\n');%say what we did here&timer

    
disp('Hello World')
end





