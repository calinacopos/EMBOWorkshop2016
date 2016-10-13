close all;
clear

%% Solving Poisson equation: Lu = u_xx + u_yy = -2(2y^3 - 3y^2 + 1) + 6(1-x^2)*(2y-1) 
%% with boundary condition u_y(x,0) = u_y(x,1) = 0
%% and u_x(0,y) = 0, u_x(1,y) = -2(2y^3 - 3y^2 + 1).
%%
%% Analytical solution: u(x,y) = (1-x^2)(2y^3-3y^2+1)

nx = 50;     
dx = 1/nx;

% Form (x,y)-space
x=linspace(0,1,nx);
y=linspace(0,1,nx);
[x_,y_]=meshgrid(x,y);

% RHS of Poisson equation
f_ = -2*(2*(y_).^3-3*(y_).^2+1)+6*(1-x_.*x_).*(2*y_-1);

e=ones(nx,1);
Lx=spdiags([e -2*e e],[-1 0 1],nx,nx);
Ly=spdiags([e -2*e e],[-1 0 1],nx,nx);

% Homoegeneous Neumann BC on y-direction
Ly(1,1) = -2;
Ly(1,2) = 2;
Ly(nx,nx-1) = 2;
Ly(nx,nx) = -2;

% Nonhomog. Neumann BC on x-direction
Lx(1,1) = -2;
Lx(1,2) = 2;
Lx(nx,nx-1) = 2;
Lx(nx,nx) = -2;
f_(:,nx) = f_(:,nx) - reshape(2*dx*(-2*(2*(y).^3-3*(y).^2+1))/(dx^2),nx,1);

% make 1D identities
Ix = speye(nx);
Iy = speye(nx);

% form 2D matrix from kron
Adiff = (kron(Iy,Lx) + kron(Ly,Ix))/(dx^2);

u_num = Adiff\reshape(f_,nx^2,1);

figure(1);
u_exact = (1-(x_.*x_)).*(2*(y_).^3-3*(y_).^2+1);
pcolor(x_,y_,u_exact); colorbar;
title('Analytical solution');

figure(2);
pcolor(x_,y_,reshape(u_num,nx,nx)); colorbar;
title('Numerical solution');