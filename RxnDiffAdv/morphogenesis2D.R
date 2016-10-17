#
# simulate logistic equation with a threshold and recovery variable in space
#   the equations are essentially the FitzHugh-Nagumo equations to
#   describe an excitable medium
#

# set parameters
# --------------

dt	= 0.001   # time step
n 	= 100   # number of space steps
T 	= 100000  # number of time steps
tplot = 1000   # output time steps

Da 	= 10   # diffusion coefficient for activator
pa 	= dt*Da  # fraction moving on each step

Dh 	= 200  #diffusion coefficient for inhibitor
ph 	= dt*Dh # fraction moving on each step
p1 = 1.0 	#production of activator
p2 = 1.0 	#production of inhibitor
m = 1.0  	#decay rate for activator (m < v)
v = 1.2  	# decay rate for inhibitor

K = 0.03;

# function to pause program
# -------------------------
pauseit <- function(x) { k <- proc.time()
Sys.sleep(x)
proc.time() - k # The cpu usage should be negligible
}


# growth function
#
fa=function(a,h) dt*(p1*(a*a/(h*(1.0+K*a*a)))-a*m)
fh=function(a,h) dt*(p2*a*a-v*h)



# set initial conditions
# ----------------------
#a   =matrix(0.5,n,n)
anew=matrix(0,n,n)
a = matrix(runif(n*n, .1, 1), n, n)
#N[1,1:n]=N0
#a[1:10,1:10]=N0

h	= matrix(0.1,n,n)
hnew = matrix(0,n,n)

# set random initial conditions
# ----------------------
#a=runif(n,min=0.1,max=1)
#h=runif(n,min=0.1,max=1)

# initilize matrix corresponding to the diffusion operator -- no flux conditions
#
Adiff = matrix(0,n,n)
Adiff[1,1]   = -pa
Adiff[1,2]   =  pa
Adiff[n,n-1] =  pa
Adiff[n,n]   = -pa

for (j in 2:(n-1) ){
  
  Adiff[j,j-1] =	pa
  Adiff[j,j  ] = -2*pa
  Adiff[j,j+1] =	pa
}

hdiff = matrix(0,n,n)
hdiff[1,1]   = -ph
hdiff[1,2]   =  ph
hdiff[n,n-1] =  ph
hdiff[n,n]   = -ph

for (j in 2:(n-1) ){
  
  hdiff[j,j-1] =	ph
  hdiff[j,j  ] = -2*ph
  hdiff[j,j+1] =	ph
}



# run simulation
# --------------
#
for (t in 1:T){
  
  #   if(t==5000){
  # 	N[40:50,40:50]=N[40:50,40:50] + N0
  #   }
  
  
  anew = a + fa(a,h) + (Adiff %*% a) + (a %*% Adiff)
  hnew = h + fh(a,h) + (hdiff %*% h) + (h %*% hdiff)
  
  # replace N with the new valuesa
  #
  a=anew
  h=hnew
  
  # make a plot of the solution
  #
  if(t %% tplot == 0) {
    #image(a,zlim=c(0,1),col=topo.colors(100),main=sprintf("time=%g",t*dt))
    image(a,zlim=c(0,1),col=grey((0:11/11)),main=sprintf("time=%g",t*dt),key=TRUE)
    #persp(N,zlim=c(-0.5,1.2))
    
    # pause
    #
    pauseit(0.001)
    
  }
  
}  # end loop in time