function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)
%optimizes the reactor


%% PRovisorisch

%idealreal 0 ideal 1 real(Peng robinson)
idealreal=0;

% INIT SOLVER
%set parameters
Pressure=101325;
FeedCH4= 0.0181;
uberschuss=1.05;
Tfeed=700;
Touter=1600;
pfrseries=2;

%idealreal=0;
% Pressure=101325;
% FeedCH4= 0.0181;
% uberschuss=1.05;
% Tfeed=700;
% Touter=1600;
% pfrseries=2  or 1;


%% ASSIGN
%Calculate and assign reactor data for one tube
rad=unt(1).rad;%Radius
l=unt(1).h;
unt(1).Aq=pi.*rad.^2;%Querschnittsfläche
unt(1).As=2.*pi.*rad.*l;%Oberfläche
unt(1).V=l.*unt(1).Aq;%Volumen
unt(1).a=unt(1).As./unt(1).V;   %specific surface (m2)/m3





%% HANDlES
%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)...
    kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)...
    PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR
cphand=@(T,cmp,unt,n) heat_capacity(T,cmp,unt,n); %handle for cp as fun of t for nth component in struct
Uhand=@(cmp,unt,p,T,F,cp,Z) HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z);%local heat transfer coefficient as function of..


%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE


MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cphand,Uhand,cmp,unt,str,idealreal);
disp('MBEB handles set')


%% SOLVE





%Integration Vector (Volume of one Tube)
Vspan=linspace(0,unt(1).V*pfrseries,100);

%Starting Values/Options
y0=[Pressure; 0; FeedCH4; uberschuss*FeedCH4; 0; 0; Tfeed; Touter];
options = odeset('NonNegative',1);

%solve
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options); %get solution


%% Assign Outputs



y=A(end,2:6)./sum(A(end,2:6));
str(5).yN2=y(1);
str(5).yCH4=y(2);
str(5).yNH3=y(3);
str(5).yH2=y(4);
str(5).yHCN=y(5);

HCNout=A(end,6);
NTubes=12.86/HCNout; %NR TUBES



%% Optimize

   function opti=optiyield(FCH4opt,MBEBhandle,options,Vspan,Pressure)       

        %Starting Values/Options
        y0opti=[Pressure; 0; FCH4opt; 1.05*FCH4opt; 0; 0; 700; 1600];
        

        [Vopt,Aopt] = ode15s(MBEBhandle,Vspan,y0opti,options);

   opti=1-Aopt(end,6)/Aopt(1,3); %find optimal yield i terms of methane, 
   %opti=1-Aopt(end,6)/sum(Aopt(end,2:6)); %find max molar fraction HCN
   end

 targetfun = @(FCH4opt)optiyield(FCH4opt, MBEBhandle,options,Vspan,Pressure);
options = optimoptions('lsqnonlin', 'display','off');
CH4_guess = FeedCH4;
CH4_min = 1e-4;
CH4_max = 0.1;
F_CH4_opti = lsqnonlin(targetfun,CH4_guess,CH4_min,CH4_max,options);

%% EVAL (to be externalized)
figure
subplot(2,2,1)
plot(Vspan,A(:,2:6))
title('F')
subplot(2,2,2)
plot(Vspan,A(:,1))
title('P')
subplot(2,2,3)
plot(Vspan,A(:,7))
title('Trxn')
subplot(2,2,4)
plot(Vspan,A(:,8))
title('Tflu')





disp('Reactor optimizer completed normally, calculated as ideal or real')
%NRTubes=12.86/HCNout



end