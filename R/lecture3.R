## -------------------------------------------- ##
## EDSD 2022-2023: Population projections
## Lecture 3
## Matrix projections
## Date: 03/05/2023
## Instructor: Ugofilippo Basellini
## -------------------------------------------- ##


## cleaning the workspace
rm(list=ls(all=TRUE))

## set working directory in active folder (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## loading useful packages
library(tidyverse)

## loading the data (derived from Preston et al. 2001)
load("data/EDSD.lecture2.Rdata")

## output from yesterday
head(my.DF)
NFx <- my.DF$NFx
sFx <- my.DF$sFx
bFx <- my.DF$bFx

## dimension of problem
m <- length(NFx)

## adjust the sFx and bFx vectors
sFx <- sFx[-m]
bFx[is.na(bFx)] <- 0
# bFx[m] <- 0

## create an empty Leslie matrix
L <- matrix(0,nrow=m,ncol=m)

## insert Sfx and bfx in Leslie matrix
L[1,] <- bFx
diag(L[-1,]) <- sFx
L[m,m] <- sFx[m-1]

## compute the projected population
NFx5 <- c(L%*%NFx)

## check that this is equal to what we got yesterday
NFx5
my.DF$NFx5
all.equal(NFx5,my.DF$NFx5)

## create a function for projecting population
## using matrix formulas
PopProj <- function(Age,AgeGroup,bFx,sFx,N0,n){
  ## dimension of problem
  m <- length(N0)
  ## adjust the sFx and bFx vectors
  sFx <- sFx[-m]
  bFx[is.na(bFx)] <- 0
  ## create the Leslie matrix
  L <- matrix(0,nrow=m,ncol=m)
  L[1,] <- bFx
  diag(L[-1,]) <- sFx
  L[m,m] <- sFx[m-1]
  ## create my output matrix
  N <- matrix(0,m,n+1)
  N[,1] <- N0
  i <- 1
  for (i in 1:n){
    N[,i+1] <- L%*%N[,i]
  }
  out <- cbind(data.frame(Age=Age,AgeGroup=AgeGroup),N)
  return(out)
}

## project 20 periods ahead
n <- 20
my.proj <- PopProj(Age=my.DF$Age,AgeGroup=my.DF$AgeGroup,bFx=my.DF$bFx,
                   sFx=my.DF$sFx,N0=my.DF$NFx,n=n)

all.equal(my.proj$`1`,my.DF$NFx)
all.equal(my.proj$`2`,NFx5)
my.proj$`21`

## long data type
my.proj.long <- my.proj %>% 
  pivot_longer(-c(Age,AgeGroup),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
    year= (period-1)*5 + 1993,
    yearF = factor(year))

## pyramid
my.proj.long %>% 
  ggplot(aes(x=AgeGroup,y=population,fill=yearF)) +
  geom_bar(data = subset(my.proj.long, year %in% c(1993,1998,2093)),
    stat="identity",position = "dodge",color="black") +
  coord_flip() +
  ggtitle("Swedish female population")






