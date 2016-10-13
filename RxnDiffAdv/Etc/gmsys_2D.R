#
# Reaction--Diffusion system of Gierer and Meinhardt taken from
#  Edelstien-Keshet book page 531.
#
#                            2
#  dx   (  x       )        d x
#  -- = ( --- - mu ) x + Dx --- 
#  dt   (  y       )          2
#                           dx
#
#
#           2                2
#  dy   (  x       )        d y
#  -- = ( --- - nu ) y + Di --- 
#  dt   (  y       )          2
#                           dx
#
#
#  homogeneous equilibrium
#        nu           nu
#   x = -----,   y = -----
#                       2
#        mu           mu
#
#  is stable when  nu > mu, without diffusion
#
#
# NOTE: I ADDED SATURATION TO THE ACTIVATOR AND THIS IS 2D
#
##########################################################################


# parameters
#
mu = 1.0     # self activation rate
nu = 1.1     # self inhibition rate
ks = 0.2     # saturation constant for activator
ks = 0.0     # saturation constant for activator
Dx = 2.0    # diffusion coefficent of activator
Dy = 20.0   # diffusion coefficent of inhibitor

#mu = 5
#nu = 6
#ks = 0.2



Nt = 50000   # number of time steps to take
dt = 0.01   # time step
tplot = 200  # make a plot every this many steps
tpause = 0.1 # system pause in animation

Nx = 50      # number of space steps 
dx = 1.0     # length of space step

# initial conditions -- equilibrium, then perturb
#
x0 = nu/mu    
y0 = nu/mu^2
X = matrix(x0 , Nx, Nx) 
Y = matrix(y0 , Nx, Nx)

# perturbation of equilibrium
#
X = X +  matrix(runif(Nx*Nx,-0.01,0.01),Nx,Nx) # perturb

# make reaction functions 
#
#Rx = function(x,y) (x  /y - mu)*x
Rx = function(x,y) (x  /(y*(1+ks*ks*x*x)) - mu)*x
Ry = function(x,y) (x^2/y - nu)*y

# per time step hoping rates from diffusion
#
px = dt*Dx/dx^2
py = dt*Dy/dx^2

# initilize matrix corresponding to the diffusion operator -- no flux conditions
#
Adiff = matrix(0,Nx,Nx)
#Adiff[1,1]     = -1.0
#Adiff[1,2]     =  1.0
#Adiff[Nx,Nx-1] =  1.0
#Adiff[Nx,Nx]   = -1.0

# periodic boundaries
#
Adiff[1,1]     = -2.0
Adiff[1,2]     =  1.0
Adiff[1,Nx]    =  1.0
Adiff[Nx,Nx-1] =  1.0
Adiff[Nx,Nx]   = -2.0
Adiff[Nx,1 ]   =  1.0

for (j in 2:(Nx-1) ){
  Adiff[j,j-1] =  1.0
  Adiff[j,j  ] = -2.0
  Adiff[j,j+1] =  1.0
}



# function to pause program 
# -------------------------
pauseit <- function(x) { p1 <- proc.time() 
                        Sys.sleep(x) 
                        proc.time() - p1 # The cpu usage should be negligible 
                       }


# begin main loop in time
#
for( t in 1:Nt ){

    Xnew = X + dt*Rx(X,Y) +  px*(Adiff %*% X) + px*(X %*% Adiff)
    Ynew = Y + dt*Ry(X,Y) +  py*(Adiff %*% Y) + py*(Y %*% Adiff)

    # update
    #
    X = Xnew
    Y = Ynew

    if(t %% tplot == 0) {

      # animate
      #
      zmax = max(X,1.005*x0)
      zmin = min(X,x0/1.005)
      filled.contour(X,zlim=c(zmin,zmax),color.palette=heat.colors)
      title(sprintf("time=%g",t*dt))
      
      # pause
      #
      pauseit(tpause)
    }
    
}  # end time loop

