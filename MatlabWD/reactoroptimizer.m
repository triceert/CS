function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)


Volmax=0.000353*800;

%initialize solver
Vspan=linspace(0,Volmax,100);
y0=[800000; 0; 0.00005; 0.00005; 0; 0; 700; 1600];


%declare needed handles
kinhand = @(T,F)kinetics(T,F); %gives rate of nth rxn
parthand  = @(p,T,F,cmp,unt,n)PRpartials(p,T,F,cmp,unt,n);%gives partial pressure for nth comp in F calc with PR

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,parthand,cmp,unt,str);
disp('MBEB handles set')

%SOLVE
options = odeset('NonNegative',1);
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0,options);


%EVAL
plot(Vspan,A(:,2:6))



disp('reactor optimizer completed normally')
end