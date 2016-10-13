% 
% simulate the FitzHugh-Nagumo equations:
% v(t) is the membrane potential (measured in mV)
% w(t) is the recovery potential (state of potassium channels)
%
%   dv        v^3
%  ---- = v - --- - w + I
%   dt         3
%
%   dw    
%  ---- = e*(v + a + b*w)
%   dt
%
% Coupled nonlinear ordinary differential equations solved with ode45.
%
% Calina Copos (Sept 30, 2016)
%

close all;
clear;

% parameters for the oscillator
% --------------
a = 0.7;
b = -0.8;
e = 0.08;
params = [a;b;e];
I = 0.5; %0.01 (small), 1 (large), 0.26 (intermediate), 500 (very large)

N = 10000; % number of time steps

% initial conditions, y0=[v0;w0]
% --------------
y0=[0.1;0];

%ode45 is a numerical integrating package that uses a variable step size
%third or fourth-order Runge-Kutta method to approximate the solution of
%the first-order ODE specified in the function 'fhn,' in this case.
%You can learn more about it by typing "help ode45" in Matlab's command
%window.
% --------------
[t,y] = ode45('fhn',linspace(0,50,N),y0,[],I,params);

% the variable y, returned by ode45, gives v(t) and w(t) as column vectors
% --------------
v  = y(:,1);
w  = y(:,2);

% plot the solutions to the coupled nonlinear ODEs
% --------------
h=figure(1); hold on; grid on; box on;
scatter(t,v,'ok','fillcolor','k');
scatter(t,w,'r','fillcolor','r');
xlim([0 50]); ylim([-2 2]);
title('FitzHugh-Nagumo system of ODEs','FontSize',30);
xlabel('Time (t)','FontSize',30);
set(gca,'FontSize',30);
L = legend('v(t): membrane potential (mV)',...
          'w(t): recovery variable');