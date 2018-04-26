function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)
%optimizes the reactor


%% PRovisorisch

%idealreal 0 ideal 1 real(Peng robinson)
idealreal=0;


%% ASSIGNn
%Calculate and assign reactor data for one tube
rad=unt(1).rad;%Radius
l=unt(1).h;
unt(1).Aq=pi.*rad.^2;%Querschnittsfläche
unt(1).As=2.*pi.*rad.*l;%Oberfläche
unt(1).V=l.*unt(1).Aq;%Volumen
unt(1).a=unt(1).As./unt(1).V;   %specific surface (m2)/m3



%% INIT SOLVER
Pressure=80325;
%FeedCH4=3.87*1e-2;
FeedCH4=0.0233;
uberschuss=1.05;
Tfeed=700;
Touter=1600;
pfrseries=6;


%Integration Vector (Volume of one Tube)
Vspan=linspace(0,unt(1).V*pfrseries,100);

%% HANDlES
%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)...
    kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)...
    PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR
cphand=@(T,cmp,unt,n) heat_capacity(T,cmp,unt,n); %handle for cp as fun of t for nth component in struct
Uhand=@(cmp,unt,p,T,F,cp,Z) HeatTransferCoefficient(cmp,unt,p,T,F,cp,Z);%local heat transfer coefficient as function of..


%Starting Values


y0=[Pressure; 0; FeedCH4; uberschuss*FeedCH4; 0; 0; Tfeed; Touter];



%% SOLVE





options = odeset('NonNegative',1);
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cphand,Uhand,cmp,unt,str,idealreal);
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options); %get solution
%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE

disp('MBEB handles set')

y=A(end,2:6)./sum(A(end,2:6))
str(5).yN2=y(1);
str(5).yCH4=y(2);
str(5).yNH3=y(3);
str(5).yH2=y(4);
str(5).yHCN=y(5);

HCNout=A(end,6);



%


opti=12.86/HCNout %NR TUBES
opti2=A(end,6)/A(1,3)

%% Assign Outputs



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