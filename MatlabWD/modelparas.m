optipars


%% WORKS 
%MODELLING PARAMETERS
unt(1).cocross=1;           % 0 cross 1 co 
                                    %Heating Medium Flow rate, only usefull when co-current flow
        str(4).G=0.06;               %in flow heating medium per fucking tube(ignored if cross heated)
unt(1).ideal_real=0;        %ideal gas 0, peng robinson 1
unt(1).nrow=6;             %number reactr elements in row

    

%THERMO
str(1).p=101325;   %feed pressure
str(1).T=800;       %feed temperature
str(2).T=1600;      %touter 

%FEEDS
str(1).FCH4=0.018; %absolute feed ch4 per single tube mol s-1   
str(1).ubsch=1.05;  %überschuss NH3


%% WORKS 2
%MODELLING PARAMETERS
unt(1).cocross=1;           % 0 cross 1 co 
                                    %Heating Medium Flow rate, only usefull when co-current flow
        str(4).G=0.09;               %in flow heating medium per fucking tube(ignored if cross heated)
unt(1).ideal_real=1;        %ideal gas 0, peng robinson 1
unt(1).nrow=10;             %number reactr elements in row

    

%THERMO
str(1).p=101325;   %feed pressure
str(1).T=800;       %feed temperature
str(2).T=1600;      %touter 

%FEEDS
str(1).FCH4=0.013; %absolute feed ch4 per single tube mol s-1   
str(1).ubsch=1.08;  %überschuss NH3


