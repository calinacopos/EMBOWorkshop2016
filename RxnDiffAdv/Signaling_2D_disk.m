%
% simulate a reaction-diffusion advection equation in a disk domain
% using a finite volume method
%
% a(t,x,y) denotes the active form of the signaling protein
% b(t,x,y) denotes the inactive form of the signaling protein
%
% Calina Copos (Sept 28, 2016)
%

close all;
clear;

% set parameters
% --------------

dt      = 0.0001;       % time step
Nt       = 200000;     % number of time steps 
tplot   = 1000;        % output time steps

N = 59;
dr = 1/(N+0.5);
dth = 2*pi/N;

r = dr*( (1:(N+1))'-0.5 );
th = dth*((1:(N+1))'-1);

N = N+1;

X = r*sin(th');
Y = r*cos(th');

Da = 0.002;              % diffusion coefficient for active          
Db = 0.4;               % diffusion coefficient for inactive

mu = 0.067;
nu = 1.0;
K = 1.0;

u1 = 0.1;            % advection flow speed in radial direction (assumed positive)
u2 = 0.002;               % zero advection flow speed in azimuthal direction

% initialize variables ( with the convention a_ij = a(r_i, theta_j) )
% --------------
a = zeros(N,N);
b = zeros(N,N);

%set initial conditions
% --------------
a = randn(N,N);
b = randn(N,N);

% reaction component
% --------------
Ra = (mu+a.^2/(K+a.^2))*b-nu*a;
Rb = -Ra;

% initialize matrix corresponding to the radial component of diffusion operator  (FD method) 
% --------------
e = ones(N,1);
rp = r+dr/2;
rm = r-dr/2;
partialLr = spdiags([rm./r -2*e rp./r],[-1 0 1],N-1,N);
Lr = partialLr(1:(N-1),1:(N-1));
Lr = [partialLr; sparse(zeros(1,N))];

% No flux on r-direction
Lr(N,N-1) = (N-1)/(N-0.5);
Lr(N,N) = -(N-1)/(N-0.5);

% Homogeneous Dirichlet in r-direction
%Lr(N,N-1) = (N-1)/(N-0.5); 
%Lr(N,N) = -2;

% initialize matrix corresponding to the angle component of diffusion operator (FD method) 
% --------------
Lth = spdiags([e -2*e e],[-1 0 1],N,N);
Irr = spdiags(1./(r.^2),0,N,N);

% -- periodic boundary conditions
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

% initialize matrix corresponding to the advection operator (1st, upwind)
% --------------
e = ones(N,1);
Dr = spdiags([-1*e e], [-1 0],N,N);

% -- no flux conditions
%Dr(1,1) = 0;
%Dr(1,2) = 0;
%Dr(N,N-1) = -1;
%Dr(N,N) = 1;

Dth = spdiags([-1*e e], [-1 0],N,N);

Ir = spdiags(1./r,0,N,N);

% run simulation
% --------------
runMode = 'A'; % mode B is animation of position 
for t=1:Nt
    Adiff = a*Lr/(dr^2) + Lth*a*Irr/(dth^2);    
    Bdiff = b*Lr/(dr^2) + Lth*b*Irr/(dth^2);
    
    Aadv = a*u1*Dr/dr + u2*Ir*a*Dth/dth;
    Badv = b*u1*Dr/dr + u2*Ir*a*Dth/dth;
    
    anew = a + Da*dt*Adiff + dt*Aadv + dt*Ra;
    bnew = b + Db*dt*Bdiff + dt*Badv + dt*Rb;
    
    a = anew;
    b = bnew;
    
    %keyboard;
    
    switch(runMode)
    case 'B'
        if mod(t,tplot) == 0 
            pcolor(X,Y,full(a));
            shading flat;
            colorbar; colormap('hot');
            caxis([-10 10]);
            title(sprintf('time=%g',t*dt))
            getframe(gcf);
        end
    end
    
    %keyboard;
end

pcolor(X,Y,full(a));
shading flat;
colorbar; colormap('hot');
%caxis([-10 10]);
title(sprintf('time=%g',t*dt))