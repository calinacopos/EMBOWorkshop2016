% lap2d.m
%
% form the (scaled) matrix for the 2D Laplacian for Dirichlet boundary
% conditions on a rectangular node-centered nx by ny grid
%
% input: nx -- number of grid points in x-direction (no bdy pts)
%        ny -- number of grid points in y-directio
%
% output: L2 -- (nx*ny) x (nx*ny) sparse matrix for discrete Laplacian


function L2 = lap2d( nx,ny )

% make 1D Laplacians
%
Lx = lap1d(nx);
Ly = lap1d(ny);

% Dirichlet BC on y-direction
%
%Ly(1,1) = -2;
%Ly(1,2) = 1;
%Ly(ny,ny-1) = 1;
%Ly(ny,ny) = -2;

% Dirichlet BC on x-direction
%
%Lx(1,1) = -2;
%Lx(1,2) = 1;
%Lx(nx,nx-1) = 1;
%Lx(nx,nx) = -2;

% Neumann BC on y-direction
Lx(1,1) = -2;
Lx(1,2) = (1+2*1/nx);
Lx(nx,nx-1) = (1+2*1/nx);
Lx(nx,nx) = -2;

% Neumann BC on y-direction
Ly(1,1) = -2;
Ly(1,2) = (1+2*1/nx);
Ly(nx,nx-1) = (1+2*1/nx);
Ly(nx,nx) = -2;

% make 1D identities
%
Ix = speye(nx);
Iy = speye(ny);

% form 2D matrix from kron
%
L2 = kron(Iy,Lx) + kron(Ly,Ix);

end


% lap1d
%
% form the (scaled) 1D Laplacian for Dirichlet boundary conditions
% on a node-centered grid
%
% input: n -- number of grid points (no bdy pts)
%
% output: L -- n x n sparse matrix for discrete Laplacian


function L = lap1d(n)

e=ones(n,1);
L = spdiags([e -2*e e],[-1 0 1],n,n);

end
