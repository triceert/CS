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

y0=[100000; 1e-25; 0.00001; 0.00001; 1e-25; 1e-25; 700; 1600];


%declare needed handles
kinhand = @(T,F,unt,PRNH3,PRCH4,idealreal)kinetics(T,F,unt,PRNH3,PRCH4,idealreal); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cmp,unt,str,idealreal);
disp('MBEB handles set')

%SOLVE
options = odeset('NonNegative',1);
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options);


%EVAL
plot(Vspan,A(:,2:6))



disp('reactor optimizer completed normally')
end