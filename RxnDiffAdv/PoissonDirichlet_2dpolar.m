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

N = 29;
dr = 1/(N+0.5); % r-coordinate centered away from r = 0 due to singularity
dth = 2*pi/N;

% Form r and theta coordinates
r = dr*( (1:(N+1))'-0.5 );
th = dth*((1:(N+1))'-1);
N = N+1;
X = r*sin(th');
Y = r*cos(th');

% RHS of Poisson equation
f_ = -ones(N,N);               %% (A)
%f_ = ones(N,N);                %% (B)
%f_ = 8.*X;                      %% (C)

e = ones(N,1);
rp = r+dr/2;
rm = r-dr/2;

partialLr = spdiags([rm./r -2*e rp./r],[-1 0 1],N-1,N);
Lr = partialLr(1:(N-1),1:(N-1));
Lr = [partialLr; sparse(zeros(1,N))];

% No flux on r-direction
%Lr(N,N-1) = (N-1)/(N-0.5);
%Lr(N,N) = -(N-1)/(N-0.5);

% Homogeneous Dirichlet in r-direction
Lr(N,N-1) = (N-1)/(N-0.5); 
Lr(N,N) = -2;
%f_(N,:) = f_(N,:) - 0.5*N/((N-0.5)*dr^2);                          %% (B)
%f_(N,:) = f_(N,:) - reshape(sin(th)*N/((N-0.5)*dr^2),1,N);          %% (C)

Lth = spdiags([e -2*e e],[-1 0 1],N,N);
Irr = spdiags(1./(r.^2),0,N,N);

% Periodic BC in theta-direction
Lth(1,1) = -2;
Lth(1,2) = 1;
Lth(1,N) = 1;
Lth(N,1) = 1;
Lth(N,N-1) = 1;
Lth(N,N) = -2;

% make 1D identities
Ith = speye(N);

% form 2D matrix from kron (kronecker command)
Adiff = kron(Ith,Lr)/(dr^2) + kron(Lth,Irr)/dth^2;

u_num = Adiff\reshape(f_,N^2,1);

% Plot analytical solution
figure(1);
u_exact = (1-X.^2-Y.^2)/4;                 %% (A)
%u_exact = (1+X.^2+Y.^2)/4;                 %% (B)
%u_exact = (X.^2+Y.^2).*X;                   %% (C)
pcolor(X,Y,u_exact); colorbar; caxis([-1 1]);
title('Analytical solution');

% Plot numerical solution
figure(2);
pcolor(X,Y,reshape(u_num,N,N)); colorbar; caxis([-1 1]);
title('Finite difference solution');