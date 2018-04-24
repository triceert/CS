
%% DATA OPENER
%Function Description: A function that opens up to date data into a struct
%Arguments:
%      table_string        file name of the excel table containing
%                               compound and unit operations data base
%           
%Outputs:
%      compstruct        the  according compound data struct
%      unitstruct        the  according compound data struct
%      streamstruct      the according streams struct

function [compstruct,unitstruct,streamstruct]=dataopener(table_string)
    compstruct=structmaker(table_string,1);
    unitstruct=structmaker(table_string,2);
    streamstruct=structmaker(table_string,3);
    disp('Dataopener terminated successfully')
end


%% %%%%%%%%%%%%%%%
%NESTED FUNCTIONS

%% STRUCTMAKER
%Function Description: A  helper function for data opener that makes...
%                      excel tables to a nice struct
%
%Arguments: string      string of excel table name            i
%             i         sheet number
%Output:    struct      The according data struct

function struct=structmaker(string,i)
    
    [~,~,ImpArr]=xlsread(string,i); % Import Excel as mixed array (Sheet i)           
            ImpArr(1,:)=[];                 %Delte Row And Col Headers
            ImpArr(:,1)=[];
            ImpArr(:,end:end-1)=[];    %delete last two columsn from excel
            ImpArr=rot90(ImpArr);           %flip for same order top as read from 
            ImpArr=flipud(ImpArr);          %left to right in excel file
    NewTab=array2table(ImpArr(2:end,:),...
    'VariableNames',ImpArr(1,:));    %get the table back and name fields properly
    struct=table2struct(NewTab);     %give back the final struct;
    disp('Structmaker terminated successfully')
end