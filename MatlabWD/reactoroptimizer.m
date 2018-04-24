function [cmp,unt,str]=reactoroptimizer(cmp,unt,str)



Vspan=linspace(0,1,100);
y0=[0; 1; 3; 0.1; 0.2; 0.3; 0.4; 0.5];

i=@(i,j) iwas(i,j);




MBEBhandle = @(t,A)MBEBpfr(t,A,i,cmp,unt,str);
disp('MBEB handle set')
[Vspan,A] = ode45(MBEBhandle,Vspan,y0)



plot(Vspan,A)



disp('reactor completed normally')
end