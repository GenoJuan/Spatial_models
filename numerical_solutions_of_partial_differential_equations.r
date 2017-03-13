#!/usr/bin/R

library(plotrix)

########################################################################################################################
########################################################################################################################

### Implement the Laplacian in a 1-D domain.
laplacian1D = function(x, neighbour.order = 1, dx, D = 1){
    x.left  = c(x[length(x)], x[1:(length(x)-1)])
    x.right = c(x[2:length(x)], x[1])
    
    y = D*(x.left + x.right - 2*x)/dx**2
    
    return(y)
}

### Implement the Lapaclian in a 2-D domain.
laplacian2D = function(x, dx, D = 1){
    x.up    = rbind(x[nrow(x),], x[1:(nrow(x)-1),])
    x.down  = rbind(x[2:nrow(x),], x[1,])
    x.left  = cbind(x[,ncol(x)], x[,1:(ncol(x)-1)])
    x.right = cbind(x[,2:ncol(x)], x[,1])
    
    y = D*(x.up + x.down + x.left + x.right - 4*x)/dx**2
}

### Swift-Hohenberg equation
SH = function(x, epsilon = 0.1){
    y = epsilon*x - x - x**3 
    return(y)
}

### Implement the Runge-Kutta order 4 numerical method
Runge.Kutta = function(x, dt = 0.1, epsilon = 0.1){
    k1 = SH(x, epsilon)*dt
    k2 = SH(x + 0.5*k1, epsilon)*dt
    k3 = SH(x + 0.5*k2, epsilon)*dt
    k4 = SH(x + k3, epsilon)*dt
    
    y = x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    return(y)
    
}

############################################################
############################################################
###                Swift-Hohenberg equation 1D           ###
############################################################
############################################################

##### Spacial parameters ######
L = 10*2*pi
dx = 0.5
x = seq(0, L, dx)

##### Temporal paramters ######
T = 100
dt = 0.001
times = seq(0, T, dt)

###### Model parameters ######

epsilon = 0.1
values = seq(-sqrt(epsilon), sqrt(epsilon), 0.01)

####### System initialization #######

space.init  = sample(x = values, size = length(x), replace = T)

####### Simulation ######

space.t = space.init

#pdf(file = "SH_1D.pdf")
plot(x, space.t, type = "l", xlab = "1-D Position", ylab = "Value")
for(time in 2:length(times)){
    space.prime = Runge.Kutta(space.t, dt = dt, epsilon = epsilon)
    space.prime.two = space.prime - 2*laplacian1D(space.prime, dx = dx, D = 1.0)*dt -  laplacian1D(laplacian1D(space.prime, dx = dx, D = 1.0), dx = dx, D = 1.0)*dt
    space.t = space.prime.two
    
    if(time %% 1000 == 0){
        plot(x, space.t, type = "l", xlab = "1-D Position", ylab = "Value")
        print(times[time])
    }
    
}

#dev.off()




############################################################
############################################################
###                Swift-Hohenberg equation 2D           ###
############################################################
############################################################

##### Spacial parameters ######
L = 10*2*pi
dx = 0.5
x = seq(0, L, dx)

##### Temporal paramters ######
T = 600
dt = 0.001
times = seq(0, T, dt)

###### Model parameters ######

epsilon = 0.1
values = seq(-sqrt(epsilon), sqrt(epsilon), 0.01)

####### System initialization #######

space.init = matrix(data = sample(x = values, size = length(x)**2, replace = T), nrow = length(x), ncol = length(x))
color2D.matplot(space.init, xlab = "X-Position", ylab = "Y-Position")
####### Simulation ######

space.t = space.init


#pdf(file = "SH_2D.pdf")
color2D.matplot(space.t, xlab = "X-Position", ylab = "Y-Position")
for(time in 2:length(times)){
    space.prime = Runge.Kutta(space.t, dt = dt, epsilon = epsilon)
    space.prime.two = space.prime - 2*laplacian2D(space.prime, dx = dx, D = 1.0)*dt -  laplacian2D(laplacian2D(space.prime, dx = dx, D = 1.0), dx = dx, D = 1.0)*dt
    space.t = space.prime.two
    
    if(time %% 1000 == 0){
        print(times[time])
        color2D.matplot(space.t, xlab = "X-Position", ylab = "Y-Position")
    }
    
}
#dev.off()
color2D.matplot(space.t, xlab = "X-Position", ylab = "Y-Position")