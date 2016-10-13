close all;
clear

%% Solving Poisson equation: Lu = u_xx + u_yy = 6xy(1-y)-2x^3 
%% with boundary condition u(x,0) = u(x,1) = 0
%% and u(0,y) = 0, u(1,y) = y(1-y).
%%
%% Analytical solution: u(x,y) = x^3 y (1-y)

nx = 50;     
dx = 1/nx;

% Form (x,y)-space
x=linspace(0,1,nx);
y=linspace(0,1,nx);
[x_,y_]=meshgrid(x,y);

% RHS of Poisson equation
f_ = 6.*x_.*y_.*(1-y_)-2*(x_).^3;

e=ones(nx,1);
Lx=spdiags([e -2*e e],[-1 0 1],nx,nx);
Ly=spdiags([e -2*e e],[-1 0 1],nx,nx);

% Homogeneous Dirichlet BC on y-direction
Ly(1,1) = -2;
Ly(1,2) = 1;
Ly(nx,nx-1) = 1;
Ly(nx,nx) = -2;

% Nonhomog. Dirichlet BC on x-direction
Lx(1,1) = -2;
Lx(1,2) = 1;
Lx(nx,nx-1) = 1;
Lx(nx,nx) = -2;
f_(:,nx) = f_(:,nx) - reshape((y.*(1-y))/(dx^2),nx,1);

% make 1D identities
Ix = speye(nx);
Iy = speye(nx);

% form 2D matrix from kron
Adiff = (kron(Iy,Lx) + kron(Ly,Ix))/(dx^2);

u_num = Adiff\reshape(f_,nx^2,1);

figure(1);
u_exact = ((x_).^3).*(y_).*(1-y_);
pcolor(x_,y_,u_exact); colorbar;
title('Analytical solution');

figure(2);
pcolor(x_,y_,reshape(u_num,nx,nx)); colorbar;
title('Finite difference solution');