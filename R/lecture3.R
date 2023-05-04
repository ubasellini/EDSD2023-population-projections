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
  ggtitle("Swedish female population") +
  theme_bw()+
  scale_fill_manual(name="Year",values=c("#E69F00", "#56B4E9","#1C7C54"))


## ---- ANIMATION INTERMEZZO ----

## 1) animated pyramids
library(viridis)
plots <- list()
my.cols <- cividis(n+1)
my.years <- unique(my.proj.long$year)
i <- 2
for (i in 1:(n+1)){
  gg <- ggplot(my.proj.long,aes(x=AgeGroup,y=population,fill=yearF)) +
    geom_bar(data = subset(my.proj.long, period == i),
             stat = "identity",color = "black") +
    coord_flip() +
    theme_bw() + ylim(0,max(my.proj.long$population)) +
    theme(legend.position = "none") +
    ggtitle(paste("Swedish female population, year",my.years[i])) +
    scale_fill_manual(values=my.cols[i])
  plots[[i]] <- gg
}
## saving plots in a single file
pdf("figs/myAnimFig.pdf")
invisible(lapply(plots, print))
dev.off()


## gganimate
library(gganimate)
library(gifski)
gg <- ggplot(my.proj.long,aes(x=AgeGroup,y=population,fill=yearF)) +
  geom_bar(stat = "identity",position = "dodge",color = "black") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis_d(option="cividis")

gg + transition_states(yearF) +
  ggtitle("Swedish female population, year {closest_state}")

anim_save("figs/F2.gif")


## alternative shiny app
my.proj.long %>% 
  plot_ly(
    x=.~AgeGroup,
    y=.~population,
    frame=.~yearF,
    type="bar"
  )


## ---- male population with matrix projections ----

## inputs
NMx <- my.DF$NMx
sMx <- my.DF$sMx
bMx <- my.DF$bMx

## adjust the sMx and bMx vectors
sMx <- sMx[-m]
bMx[is.na(bMx)] <- 0


## create a new Leslie matrix

L <- matrix(0,2*m,2*m)

## compute BM and LM
LM <- BM <- LF <- matrix(0,nrow=m,ncol=m)

## insert Sfx and bfx in Leslie matrix
BM[1,] <- bMx
diag(LM[-1,]) <- sMx
LM[m,m] <- sMx[m-1]

## for females
LF[1,] <- bFx
diag(LF[-1,]) <- sFx
LF[m,m] <- sFx[m-1]

## insert the LM, LF and BM into L
L[1:m,1:m] <- LF
L[1:m + m,1:m] <- BM
L[1:m + m,1:m + m] <- LM

## second version of L
ZEROS <- matrix(0,m,m)
LUP <- cbind(LF,ZEROS)
LDOWN <- cbind(BM,LM)
L2 <- rbind(LUP,LDOWN)
all.equal(L,L2)

## project the female+male population
Nx5 <- c(L%*%c(NFx,NMx))
all.equal(Nx5[1:m],my.DF$NFx5)
all.equal(Nx5[1:m+m],my.DF$NMx5)


## create a function for projecting population
## using matrix formulas
PopProj <- function(Age,AgeGroup,bFx,bMx,sFx,sMx,NF0,NM0,n){
  ## dimension of problem
  m <- length(NF0)
  ## adjust the sFx and bFx vectors
  sFx <- sFx[-m]
  bFx[is.na(bFx)] <- 0
  ## adjust the sMx and bMx vectors
  sMx <- sMx[-m]
  bMx[is.na(bMx)] <- 0
  ## create the Leslie matrix
  L <- matrix(0,nrow=2*m,ncol=2*m)
  LM <- BM <- LF <- matrix(0,nrow=m,ncol=m)
  ## insert Sfx, bfx, bmx, smx in Leslie matrix
  BM[1,] <- bMx
  diag(LM[-1,]) <- sMx
  LM[m,m] <- sMx[m-1]
  LF[1,] <- bFx
  diag(LF[-1,]) <- sFx
  LF[m,m] <- sFx[m-1]
  
  ## insert the LM, LF and BM into L
  L[1:m,1:m] <- LF
  L[1:m + m,1:m] <- BM
  L[1:m + m,1:m + m] <- LM
  
  ## create my output matrix
  N <- matrix(0,2*m,n+1)
  N[,1] <- c(NF0,NM0)
  i <- 1
  for (i in 1:n){
    N[,i+1] <- L%*%N[,i]
  }
  out <- cbind(data.frame(Age=Age,AgeGroup=AgeGroup,
                          Sex=c(rep("Females",m),rep("Males",m))),N)
  return(out)
}


## project 20 periods ahead
n <- 20
my.proj <- PopProj(Age=my.DF$Age,AgeGroup=my.DF$AgeGroup,bFx=my.DF$bFx,
                   sFx=my.DF$sFx,bMx=my.DF$bMx,
                   sMx=my.DF$sMx,NF0=my.DF$NFx,NM0=my.DF$NMx,n=n)

all.equal(my.proj$`1`,c(my.DF$NFx,my.DF$NMx))
all.equal(my.proj$`2`,c(my.DF$NFx5,my.DF$NMx5))

## long data type
my.proj.long <- my.proj %>% 
  pivot_longer(-c(Age,AgeGroup,Sex),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
         year= (period-1)*5 + 1993,
         yearF = factor(year))

## pyramid
my.proj.long %>% 
  ggplot(aes(x=AgeGroup,y=population,fill=yearF)) +
  geom_bar(data = subset(my.proj.long, year %in% c(1993,1998,2093) & Sex == "Males"),
           stat="identity",position = "dodge",color="black",mapping = aes(y = -population)) +
  geom_bar(data = subset(my.proj.long, year %in% c(1993,1998,2093) & Sex == "Females"),
           stat="identity",position = "dodge",color="black") +
  coord_flip() +
  ggtitle("Swedish female population") +
  theme_bw()+
  scale_fill_manual(name="Year",values=c("#E69F00", "#56B4E9","#1C7C54")) +
  scale_y_continuous(limits=c(-3.5e5,3.5e5),
                     breaks = seq(-4e5,4e5,1e5),
                     labels = abs(seq(-4e5,4e5,1e5))) +
  geom_text(data = subset(my.proj.long, period %in% c(1)),
            aes(y = max(population)/1.25, x = 17, label="Females"),size=7) +
  geom_text(data = subset(my.proj.long, period %in% c(1)),
            aes(y = -max(population)/1.25, x = 17, label="Males"),size=7)



## saving the data for tomorrow
save.image("data/EDSD.lecture3.Rdata")

## END