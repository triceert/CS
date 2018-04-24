function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)


%initialize solver
Vspan=linspace(0,1,100);
y0=[1; 1; 0; 0.1; 0.2; 0.3; 0.4; 0.5];


%declare needed handles
kinhand = @(T,pCH4,pNH3,n)kinetics(T,pCH4,pNH3,n);
prhand  = @()

%MAIN HANDLE CONTAINING ALL OTHER HANDLES FROM ABOVE
MBEBhandle = @(t,A)MBEBpfr(t,A,kinhand,prhand,cmp,unt,str);
disp('MBEB handles set')

%SOLVE
[Vspan,A] = ode15s(MBEBhandle,Vspan,y0);


%EVAL
plot(Vspan,A)



disp('reactor optimizer completed normally')
end