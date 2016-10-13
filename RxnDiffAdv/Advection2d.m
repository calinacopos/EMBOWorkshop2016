close all;
clear

%% Solving u_t + v.(u_x + u_y) = 0 
%% with initial condition u(x,y,0) = exp(-((x-0.5)^2+y^2)/L^2)

dt = 0.001;     % time step
Nt = 10000;     % number of time steps
nx = 50;     
dx = 1/nx;
v = 1; w = 0.01; L = 2;

dv = v*dt/dx;
dw = w*dt/dx;

% Form (x,y)-space
x=linspace(0,1,nx);
y=linspace(0,1,nx);
[x_,y_]=meshgrid(x,y);

% RHS of Poisson equation
f_ = exp(-((x_-0.5).^2+(y_).^2)/(L^2));

e=ones(nx,1);
Dx=spdiags([-1*e e],[-1 0],nx,nx);
Dy=spdiags([-1*e e],[-1 0],nx,nx);

u = zeros(nx,nx);
u = f_; % initial condition

pcolor(full(f_));
colorbar; caxis([0.5 1]);
title('Initial condition, T=0');
pause;

for t=1:Nt
    unew = u - dv*u*Dx - dw*Dy*u;
    u = unew;
end

pcolor(full(u));
colorbar; caxis([-1 1]);
title('Finite difference solution, T = 10,000');