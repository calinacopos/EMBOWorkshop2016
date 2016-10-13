close all;
clear

%% solving u_xx + u_yy = -2(2y^3 - 3y^2 + 1) + 6(1-x^2)*(2y-1) 
%% with boundary condition u_y(x,0)=u_y(x,1)=0 and u(0,y)=0, u(1,y) = (2y^3-3y^2+1)

nx = 50;     
dx = 1/nx;

x=linspace(0,1,nx);
y=linspace(0,1,nx);
[x_,y_]=meshgrid(x,y);
f_ = -2*(2*(y_).^3-3*(y_).^2+1)+6*(1-x_.*x_).*(2*y_-1);
u_anal = (1-(x_.*x_)).*(2*(y_).^3-3*(y_).^2+1);

e=ones(nx,1);
Lx=spdiags([e -2*e e],[-1 0 1],nx,nx);
Ly=spdiags([e -2*e e],[-1 0 1],nx,nx);

% Homogeneous Neumann BC on y-direction
Ly(1,1) = -2;
Ly(1,2) = 2;
Ly(nx,nx-1) = 2;
Ly(nx,nx) = -2;

% Nonhomog. Dirichlet BC on x-direction
Lx(1,1) = -2;
Lx(1,2) = 1;
Lx(nx,nx-1) = 1;
Lx(nx,nx) = -2;
f_(:,nx) = f_(:,nx) - reshape((2*y.^3-3*y.^2+1)/(dx^2),nx,1);

% make 1D identities
Ix = speye(nx);
Iy = speye(nx);

% form 2D matrix from kron
Adiff = (kron(Iy,Lx) + kron(Ly,Ix))/(dx^2);

u_num = Adiff\reshape(f_,nx^2,1);

figure(1);
pcolor(x_,y_,reshape(u_num,nx,nx)); colorbar;

figure(2);
pcolor(x_,y_,u_anal); colorbar;