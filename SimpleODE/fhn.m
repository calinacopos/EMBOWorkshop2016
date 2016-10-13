function dV = fhn(t,V,options,I,params)
%This is a function that goes into one of Matlab's ODE solvers.  In
%general, it takes a vector V and a time t as input and returns the time
%derivative of V. Here we are solving the FitzHugh-Nagumo system of 
%ODEs where v is the membrane potential and w is the recovery variable
%(state of potassium channels). We can use one of Matlab's solvers to
%solve these ODEs numerically.
a = params(1);
b = params(2);
e = params(3);

v = V(1); 
w = V(2);

dv = v-(v^3)/3-w+I;
dw = e*(v+a+b*w); 

dV=[dv; dw];
