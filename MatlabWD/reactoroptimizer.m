function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)
%optimizes the reactor


%% INITIAL VALUES

%idealreal 0 ideal 1 real(Peng robinson)



%set parameters BEST FOR REAL
% idealreal=1;
% Pressure=101325; %pressure
% FeedCH4= 0.0022; %Feed Methan
% uberschuss=1.05; %Feedfactor NH3
% Tfeed=700;       %Inlet Feed Temperature
% Touter=1600;     %Rxn Mixture Temperature
% pfrseries=2;     %how many PFRs in Series

% idealreal=0; %best for ideal
% Pressure=101325;
% FeedCH4=  0.0590;
% uberschuss=1.05;
% Tfeed=700;
% Touter=1600;
% pfrseries=  1;


idealreal=0; %best for ideal
Pressure=101325;
FeedCH4= 0.0181;
uberschuss=1.05;
Tfeed=700;
Touter=1600;
pfrseries=  1;


%% ASSIGN
% assign or get reactor data for one tube
rad=unt(1).rad;%Radius
l=unt(1).h;%Länge
unt(1).Aq=pi.*rad.^2;%Querschnittsfläche
unt(1).As=2.*pi.*rad.*l;%Oberfläche
unt(1).V=l.*unt(1).Aq;%Volumen
unt(1).a=unt(1).As./unt(1).V;   %specific surface (m2)/m3
MW_in = extractfield(cmp(2:6),'MW')';
R=unt(5).idgc;
HCNneeded=12.86;

%calc additional  reactor data dependent from IC and assign
FeedNH3=uberschuss*FeedCH4;
Ftot_in=FeedCH4+FeedNH3;
MW_mix_in =  (FeedCH4.*MW_in(2)+FeedNH3*MW_in(3))./Ftot_in;    
yCH4=FeedCH4/Ftot_in;
yNH3=FeedNH3/Ftot_in;
[~,Zin] = PRpartials(Pressure,Tfeed,[0;FeedCH4;FeedNH3;0;0],cmp,unt,1);
Q_in = Ftot_in*Zin*R*Tfeed/Pressure;
rho_mix_in = Zin*R*Tfeed/(Pressure*MW_mix_in); %[m3.kg-1]


str(1).p=Pressure;
str(1).T=Tfeed;
str(2).T=Touter;
str(1).yCH4=yCH4;
str(1).yNH3=yNH3;
unt(1).Q_in=Q_in;
unt(1).rho_mix_in=rho_mix_in;



%% HANDlES
%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)...
    kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)...
    PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR
cphand=@(T,cmp,unt,n) heat_capacity(T,cmp,unt,n); %handle for cp as fun of t for nth component in struct
Uhand=@(cmp,unt,p,T,F,cp,Z,heatedwith) HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z,heatedwith);%local heat transfer coefficient as well as flow and Reynolds number inside reactor




%% SOLVE

%Integration Vector (Volume of one Tube)
Vspan=linspace(0,unt(1).V*pfrseries,100);

%Starting Values/Options
y0=[Pressure; 0; FeedCH4; uberschuss*FeedCH4; 0; 0; Tfeed; Touter;0];
options = odeset('NonNegative',1);

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cphand,Uhand,cmp,unt,str,idealreal);
disp('MBEB handles set')


%solve
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options); %get solution


%% Assign Outputs



y=A(end,2:6)./sum(A(end,2:6));
str(5).yN2=y(1);
str(5).yCH4=y(2);
str(5).yNH3=y(3);
str(5).yH2=y(4);
str(5).yHCN=y(5);



Qneeded=A(end,9);


HCNout=A(end,6);
NTubes=HCNneeded/HCNout; %NR TUBES
unt(1).N_tubes=NTubes.*pfrseries;
str(5).G=sum(A(end,2:6))*NTubes;

%Assign Stream corrected with n tubes
str(1).G=Ftot_in*NTubes;
unt(1).Q_tot=Qneeded.*NTubes;


%% Optimize

   function opti=optiyield(FCH4opt,MBEBhandle,options,Vspan,Pressure)       

        %Starting Values/Options
        y0opti=[Pressure; 0; FCH4opt; 1.05*FCH4opt; 0; 0; 700; 1600;0];
        

        [~,Aopt] = ode15s(MBEBhandle,Vspan,y0opti,options);
        
        NTubesopt=12.86/Aopt(end,6); %NR TUBES
    

   %opti=1-Aopt(end,6)/Aopt(1,3); %find optimal yield i terms of methane, 
   opti=NTubesopt./(Aopt(end,6)/Aopt(1,3));
  
   end

 targetfun = @(FCH4opt)optiyield(FCH4opt, MBEBhandle,options,Vspan,Pressure);
options = optimoptions('lsqnonlin', 'display','off');
CH4_guess = FeedCH4;
CH4_min = 1e-4;
CH4_max = 0.1;
F_CH4_opti = lsqnonlin(targetfun,CH4_guess,CH4_min,CH4_max,options);

%% EVAL (to be externalized)
figure
subplot(4,1,1)
plot(Vspan,A(:,2:6))
title('F')
subplot(4,1,2)
plot(Vspan,A(:,1))
title('P')
subplot(4,1,3)
plot(Vspan,A(:,7))
title('Trxn')
subplot(4,1,4)
plot(Vspan,A(:,8))
title('Tflu')




disp('Reactor optimizer completed normally, calculated as ideal or real')




end