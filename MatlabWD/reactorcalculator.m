function [cmp,unt,str]=reactorcalculator(cmp,unt,str,plotpar)
%Function for Reactor Calculaion and Manual optimization if wanted
%% GET INITIAL VALUES






%normal inputs defined bevore

Pressure=str(1).p;
Tfeed=str(1).T;
Touter=str(2).T;
FeedCH4=str(1).FCH4;%absolute feed ch4 per single tube mol s-1    
idealreal=unt(1).ideal_real;%
uberschuss=str(1).ubsch;
pfrseries=unt(1).nrow;


%% ASSIGN FOR ONE TUBE
% assign or get reactor data for one tube
rad=unt(1).rad;%Radius
l=unt(1).h;%Länge
unt(1).Aq=pi.*rad.^2;%Querschnittsfläche
unt(1).As=2.*pi.*rad.*l;%Oberfläche
unt(1).V=l.*unt(1).Aq;%Volumen
unt(1).a=unt(1).As./unt(1).V;   %specific surface (m2)/m3
MW_in = extractfield(cmp(2:6),'MW')';
R=unt(5).idgc;
HCNneeded=unt(1).HCNneeded;

%calc additional  reactor data dependent from IC and assign
%inlet Feed calculations 
FeedNH3=uberschuss*FeedCH4;
Ftot_in=FeedCH4+FeedNH3;
MW_mix_in =  (FeedCH4.*MW_in(2)+FeedNH3*MW_in(3))./Ftot_in;    
yCH4=FeedCH4/Ftot_in;
yNH3=FeedNH3/Ftot_in;
str(1).yCH4=yCH4;
str(1).yNH3=yNH3;

R=unt(5).idgc; %ideal gas constant

%thermo calculation and assignment for IC

            switch idealreal
                case 1 %use Peng Robinson
                [~,Zin] = PRpartials(Pressure,Tfeed,[0;FeedCH4;FeedNH3;0;0],cmp,unt,1);


                case 0
                    Zin=1;

                otherwise
                    error('Wrong idealreal paramter input, use 0 or 1')

            end
        
    Q_in = Ftot_in*Zin*R*Tfeed/Pressure; %volumetric flow for EBMB
    unt(1).Q_in=Q_in;%assign volumetric flow for EBMB


%mixture rho of inlet
rho_mix_in = Zin*R*Tfeed./(Pressure*MW_mix_in); %[m3.kg-1]
unt(1).rho_mix_in=rho_mix_in;



%% HANDlES%INTEGRATION LIMITS
%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)...
    kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)...
    PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR
cphand=@(T,cmp,unt,n) heat_capacity(T,cmp,unt,n); %handle for cp as fun of t for nth component in struct
Uhand=@(cmp,unt,p,T,F,cp,Z,heatedwith) HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z,heatedwith);%local heat transfer coefficient as well as flow and Reynolds number inside reactor


%Integration Vector (Volume of one Tube)
numsolrows=100; %assign how many calculation points
Vspan=linspace(0,unt(1).V*pfrseries,numsolrows);


%% SOLVE

%Starting Values/Options
y0=[Pressure; 0; FeedCH4; uberschuss*FeedCH4; 0; 0; Tfeed; Touter;0];
options = odeset('NonNegative',1);

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cphand,Uhand,cmp,unt,str,idealreal);


%solve
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options); %get solution


%% MANUAL OPTIMIZER



%                         %% Manual Optimizer
% 
%                            function opti=optiyield(FCH4opt,MBEBhandle,options,Vspan,Pressure)       
% 
%                                 %Starting Values/Options
%                                 y0opti=[Pressure; 0; FCH4opt; 1.05*FCH4opt; 0; 0; 700; 1600;0];
% 
% 
%                                 [~,Aopt] = ode15s(MBEBhandle,Vspan,y0opti,options);
% 
%                                 NTubesopt=12.86/Aopt(end,6); %NR TUBES
% 
% 
%                            opti=1-Aopt(end,6)/Aopt(1,3); %find optimal yield i terms of methane, 
%                            %opti=NTubesopt./(Aopt(end,6)/Aopt(1,3)).^5;
%                           %opti=-NTubesopt;
%                            end
% 
%                         targetfun = @(FCH4opt)optiyield(FCH4opt, MBEBhandle,options,Vspan,Pressure);
%                         options = optimoptions('lsqnonlin', 'display','off');
%                         CH4_guess = FeedCH4;
%                         CH4_min = 1e-3;
%                         CH4_max = 0.1;
%                         F_CH4_opti = lsqnonlin(targetfun,CH4_guess,CH4_min,CH4_max,options)

%% EVAL (to be externalized)
%%ASSIGN DATA

%Molar fraction assignment
y=A(end,2:6)./sum(A(end,2:6));
str(5).yN2=y(1);
str(5).yCH4=y(2);
str(5).yNH3=y(3);
str(5).yH2=y(4);
str(5).yHCN=y(5);

%yield in terms of CH4
yield=A(end,6)/A(1,3);
unt(1).yield=yield;




%tubes calculation and assignment
HCNout=A(end,6);
NTubesaside=HCNneeded/HCNout; %NR TUBES
unt(1).N_tubes_aside=NTubesaside; %number of tubes aside
unt(1).N_tubes=NTubesaside.*pfrseries; %number of tubes total


%TOTAL STREAM ASSIGNMENT
%Assign Stream corrected with n tubes
str(1).G=Ftot_in*NTubesaside; %for tubes aside Strem in
str(5).G=sum(A(end,2:6))*NTubesaside; %for tubes aside stream out

%Heat needed
Qneeded=A(end,9);
unt(1).Q_tot=Qneeded.*NTubesaside;

%conversion CH4
unt(1).conv=(str(1).G*str(1).yCH4-str(5).G*str(5).yCH4)/(str(1).G*str(1).yCH4);


%Velocity calculations assignment (FOR/FROM ONE TUBE)
    %get compressibility and thus molar volume at every point
       
        
      Videal=R.*A(:,7)./A(:,1);
      
      
      
      switch idealreal
          case 1
              Z=zeros(numsolrows,1);
              for i=1:numsolrows
              [~,Z(i,1)] = PRpartials(A(i,1),A(i,7),A(i,2:6),cmp,unt,1);
              end
          case 0
              Z=1;
      end
                                 
      Vreal=Z.*Videal;    
    
    
    %get total molar flow at every point
        Flowmatrix=(A(:,2:6));
        Flowsum=sum(Flowmatrix,2);
    %get volumetric flow and divide by area for spatial velcoity
        velo=(Flowsum.*Vreal)./unt(1).Aq;


%PLOT

switch plotpar
    
    case 1
    figure
    title('Optimized Paramters Run Output')
    subplot(5,1,1)
    plot(Vspan,A(:,2:6))
    legend('on')
    title('F')
    subplot(5,1,2)
    plot(Vspan,A(:,1))
    title('P')
    subplot(5,1,3)
    plot(Vspan,A(:,7))
    title('Trxn')
    subplot(5,1,4)
    plot(Vspan,A(:,8))
    title('Tflu')

    subplot(5,1,5)
    plot(Vspan,velo)
    title('Velprofile')
    
    disp('Reactor calculator completed normally for one run, plots no done')

    case 0
        
       
end








end