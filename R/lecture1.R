## -------------------------------------------- ##
## EDSD 2022-2023: Population projections
## Lecture 1
## The constant exponential growth model
## Date: 02/05/2023
## Instructor: Ugofilippo Basellini
## -------------------------------------------- ##


## cleaning the workspace
rm(list=ls(all=TRUE))

## starting populations & starting inputs
N0.a <- 50
N0.b <- 35
r.a <- 5/100
r.b <- 15/100
t <- seq(0,10,by=0.01)

## function for constant exponential growth model
## depending on crude growth rate r
PopProj <- function(N0,r,t){
  NT <- N0*exp(r*t)
  return(NT)
}

## compute the projected populations
NT.a <- PopProj(N0=N0.a,r=r.a,t=t)
plot(t,NT.a,t="l",lwd=2)
NT.b <- PopProj(N0=N0.b,r=r.b,t=t)
plot(t,NT.b,t="l",lwd=2)

t.hat <- t[which(NT.b>NT.a)[1]]

plot(t,NT.a,t="l",lwd=2,ylim=range(NT.a,NT.b),
     xlab="time (years)",ylab="population")
lines(t,NT.b,col=2,lwd=2)
abline(v=t.hat,lty=2)
legend("topleft",c("pop A","pop B"),col=1:2,lwd=2)

##---- including birth, death and migration crude rates
b.a <- 3/100
d.a <- 5/100
m.a <- 1/100
b.b <- 7/100
d.b <- 4/100
m.b <- 3/100

## function for constant exponential growth model
## depending on b, d and m
PopProj <- function(N0,b,d,m,t){
  r <- b - d + m
  NT <- N0*exp(r*t)
  return(NT)
}

## compute the projected populations
NT.a <- PopProj(N0=N0.a,b=b.a,d=d.a,m=m.a,t=t)
plot(t,NT.a,t="l",lwd=2)
NT.b <- PopProj(N0=N0.b,b=b.b,d=d.b,m=m.b,t=t)
plot(t,NT.b,t="l",lwd=2)

(t.hat <- t[which(NT.b>NT.a)[1]])

plot(t,NT.a,t="l",lwd=2,ylim=range(NT.a,NT.b),
     xlab="time (years)",ylab="population")
lines(t,NT.b,col=2,lwd=2)
abline(v=t.hat,lty=2)


##---- moving away from net migration rates
N0 <- 50
b <- 6/100
d <- 4/100
e <- 1/100
I <- 1.5
t <- 0:20
m <- length(t)

## function for constant exponential growth model
## depending on b, d, e and I
PopProj <- function(N0,b,d,e,I,t){
  m <- length(t)
  r <- b - d - e
  ## create an empty matrix to contain all projected values of NT
  NT <- rep(NA,m)
  NT[1] <- N0
  ## for loop
  i <- 2
  for (i in 2:m){
    NT[i] <- NT[i-1]*exp(r) + I
  }
  return(NT)
}

NT <- PopProj(N0=N0,b=b,d=d,e=e,I=I,t=t)
plot(t,NT,t="l",lwd=2)
NT[m]

##---- incorportaing future assumptions

## assumption about e
e2 <- c(rep(e,11),rep(e*2,10))
plot(t,e2)

## assumption about I
I2 <- c(rep(I,11),seq(I,0,length.out=10))
plot(t,I2)

## function for constant exponential growth model
## depending on b, d, e and I
## with time-specific assumptions
PopProj <- function(N0,b,d,e,I,t){
  m <- length(t)
  ## create an empty matrix to contain all projected values of NT
  NT <- rep(NA,m)
  NT[1] <- N0
  ## for loop
  i <- 2
  for (i in 2:m){
    r <- b - d - e[i]
    NT[i] <- NT[i-1]*exp(r) + I[i]
  }
  return(NT)
}

NT_scenario <- PopProj(N0=N0,b=b,d=d,e=e2,I=I2,t=t)


plot(t,NT,t="l",lwd=2,ylim=range(NT,NT_scenario),
     xlab="time (years)",ylab="population")
lines(t,NT_scenario,col=2,lwd=2)
abline(v=10,lty=2)











