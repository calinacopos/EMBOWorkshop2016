%
% simulate a reaction-diffusion advection equation in a rectangular domain
% using a finite difference method
%
% a(t,x,y) denotes the active form of the signaling protein
% b(t,x,y) denotes the inactive form of the signaling protein
%
% Calina Copos (Sept 27, 2016)
%

close all;
clear;

% set parameters
% --------------

dt      = 0.0001;       % time step
nx      = 50;          % number of space steps
Nt      = 100000;     % number of time steps 
tplot   = 1000;        % output time steps

dx = 1/nx;              % spacing

Da = 0.002;%0.01;              % diffusion coefficient for active          
Db = 0.4;%0.2;               % diffusion coefficient for inactive

mu = 0.067;
nu = 1.0;

v = 0.1;               % advection flow speed in x-direction (assumed positive)
w = 0.002;               % advection flow speed in y-direction (assumed positive)

% per time step hoping rates from diffusion
% --------------
pa = dt*Da/dx^2;
pb = dt*Db/dx^2;

% per time step hoping rates from advection
% --------------
dv = v*dt/dx;
dw = w*dt/dx;

%% Laplacian operator
% --------------
rhs=zeros(nx,nx);
e=ones(nx,1);
Lx=spdiags([e -2*e e],[-1 0 1],nx,nx);
Ly=spdiags([e -2*e e],[-1 0 1],nx,nx);

% No-flux BC on y-direction
Ly(1,1) = -2;
Ly(1,2) = 2;
Ly(nx,nx-1) = 2;
Ly(nx,nx) = -2;

% No-flux BC on x-direction
Lx(1,1) = -2;
Lx(1,2) = 2;
Lx(nx,nx-1) = 2;
Lx(nx,nx) = -2;

% make 1D identities
Ix = speye(nx);
Iy = speye(nx);

% form 2D matrix from kron
Diff = (kron(Iy,Lx) + kron(Ly,Ix));

%% Avection operator (1st, upwind)
% --------------
e = ones(nx,1);
Dx = spdiags([-1*e e], [-1 0],nx,nx);

% -- no flux conditions
%Dx(1,1) = 0;
%Dx(1,2) = 0;
%Dx(nx,nx-1) = -1;
%Dx(nx,nx) = 1;

Dy = spdiags([-1*e e], [-1 0],nx,nx);
%Dy(1,1) = 0;
%Dy(1,2) = 0;
%Dy(nx,nx-1) = -1;
%Dy(nx,nx) = 1;

% initialize variables
% --------------
a = zeros(nx,nx);
b = zeros(nx,nx);

%set initial conditions
% --------------
a = randn(nx,nx)+1;
b = 2*(randn(nx,nx)+1);

%% Reaction component
% --------------
Ra = (mu+(a.^2)./(1+a.^2)).*b-nu*a;
Rb = -Ra;

%% run simulation
% --------------
x = linspace(0,1,nx);
y = linspace(0,1,nx);
[X,Y]=meshgrid(x,y);
runMode = 'B'; % mode B is animation of position 
for t=1:Nt    
    anew = a + pa*reshape(Diff*reshape(a,nx^2,1),nx,nx) + dt*Ra - (dv*a*Dx + dw*Dy*a);
    bnew = b + pb*reshape(Diff*reshape(b,nx^2,1),nx,nx) + dt*Rb - (dw*b*Dx + dw*Dy*b);
    
    % update
    a = anew;
    b = bnew;
    
    switch(runMode)
    case 'B'
        if mod(t,tplot) == 0 
            pcolor(X,Y,full(a));
            shading flat;
            colorbar; colormap('hot');
            caxis([-1 1]);
            title(sprintf('time=%g',t*dt))
            getframe(gcf);
            %pause;
        end
    end
end

%% plot end results
% --------------
pcolor(X,Y,full(a));
shading flat;
colorbar; colormap('hot');
title(sprintf('time=%g',t*dt))
caxis([-1 1]);

