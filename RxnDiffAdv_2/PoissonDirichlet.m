close all;
clear

%% (A) Solving Poisson equation: Lu = (1/r)u_r + u_rr + (1/r^2) u_tt = -1
%% with boundary condition u(r,t) = 0 on r = 1.
%%
%% Analytical solution: u(r,t) = 0.25*(1-r^2)

%% (B) Solving Poisson equation: Lu = (1/r)u_r + u_rr + (1/r^2) u_tt = 1
%% with boundary condition u(r,t) = 1/2 on r = 1.
%%
%% Analytical solution: u(r,t) = 0.25*(1+r^2)

%% (C) Solving Poisson equation: Lu = (1/r)u_r + u_rr + (1/r^2) u_tt = 8r*sint(t)
%% with boundary condition u(r,t) = sin(t) on r = 1.
%%
%% Analytical solution: u(r,t) = r^3*sin(t)

N = 50;
dx = 2/(N-1);

% Form (x,y)-space
x = linspace(-1,1,N);
y = linspace(-1,1,N);
[X,Y] = meshgrid(x,y);
f_ = -ones(N,N);               %% (A)
%f_ = ones(N,N);                %% (B)
%f_ = 8*X;                       %% (C)

% Find the matrix of the modified discrete Laplacian
[L, g, phi] = make_matrix_rhs_circleproblem(N);
%[L, g, phi] = make_matrix_rhs_ellipseproblem(N);
g_ = reshape(g,N,N);
f_ = f_ + g_/(dx^2); % Adjust RHS

L = L/(dx^2);

u_num = L\reshape(f_,N^2,1);

% Plot analytical solution
figure(1);
u_exact = (1-X.^2-Y.^2)/4;     %% (A)
%u_exact = (1+X.^2+Y.^2)/4;     %% (B)
%u_exact = (X.^2+Y.^2).*X;       %% (C)
u_exact = u_exact.*(phi<0);
pcolor(X,Y,u_exact); colorbar; caxis([-1 1]);
title('Analytical solution');

% Plot numerical solution
figure(2);
pcolor(X,Y,reshape(u_num,N,N)); colorbar;
title('Finite difference solution');