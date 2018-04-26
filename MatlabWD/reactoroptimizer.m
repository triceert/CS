function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)

%% PRovisorisch

%idealreal 0 ideal 1 real(Peng robinson)
idealreal=0


%% ASSIGN
%Calculate and assign reactor data for one tube
unt(1).Aq=unt(1).rad.^2*pi;%Querschnittfl
unt(1).As=unt(1).h*2*pi*unt(1).rad; %Radius
unt(1).V=unt(1).Aq.*unt(1).h;%Volume
unt(1).a=unt(1).As/unt(1).V;%Specific surface of reactor



%% INIT SOLVER
%Integration Vector (Volume of one Tube)
Vspan=linspace(0,unt(1).V,100);

%Starting Values
y0=[100000; 1e-35; 0.00001; 0.00001; 1e-35; 1e-35; 700; 1600];

%% HANDlES
%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)...
    kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)...
    PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cmp,unt,str,idealreal);
disp('MBEB handles set')

%% SOLVE
options = odeset('NonNegative',1);
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options);

y=A(end,2:6)./sum(A(end,2:6))

str(5).yN2=y(1);
str(5).yCH4=y(2);
str(5).yNH3=y(3);
str(5).yH2=y(4);
str(5).yHCN=y(5);

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
end