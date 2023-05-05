## -------------------------------------------- ##
## EDSD 2022-2023: Population projections
## Lecture 4
## Matrix projections: extensions
## Date: 05/05/2023
## Instructor: Ugofilippo Basellini
## -------------------------------------------- ##


## cleaning the workspace
rm(list=ls(all=TRUE))

## set working directory in active folder (R-studio command)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load the migest package
library(migest)
library(tidyverse)

## check the fundamental RC parameters
rc_model_fund
my.pars <- rc_model_fund %>%
  select(value) %>% pull()

## function to construct RC schedule
mxRG <- function(x,pars){
  t1 <- pars[1]*exp(-pars[2]*x)
  t2 <- pars[3]*exp(-pars[4]*(x-pars[5])-exp(-pars[6]*(x-pars[5])))
  mx <- t1+t2+pars[7]
  mx <- mx/sum(mx)
  return(mx)
}
## five-year age groups (works well also with one-year)
x <- seq(0,85,5)
mx <- mxRG(x=x,pars=my.pars)
plot(x, mx, type="o",pch=16)
## assume a total of 100000 net migration counts
I <- 1e5
Ix <- I*mx
sum(Ix)
plot(x, Ix, type="o",pch=16,
     xlab = "Age group",ylab= "Net migrant counts",
     main="RC migration schedule for 100,000 net migrants")


## loading yesterday's data 
load("data/EDSD.lecture3.Rdata")

I <- 25000
Ix <- I*mx
plot(x, Ix, type="o",pch=16,
     xlab = "Age group",ylab= "Net migrant counts",
     main="RC migration schedule for 25,000 net migrants")

## include net migration in our population projection function
## create a function for projecting population
## using matrix formulas
PopProjWithMigration <- function(Age,AgeGroup,bFx,sFx,Ix,N0,n){
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
    N[,i+1] <- L%*%(N[,i]+Ix/2) + Ix/2
  }
  out <- cbind(data.frame(Age=Age,AgeGroup=AgeGroup),N)
  return(out)
}

n <- 20
my.proj.nomigr <- PopProj(Age=my.DF$Age,AgeGroup=my.DF$AgeGroup,bFx=my.DF$bFx,
                          sFx=my.DF$sFx,N0=my.DF$NFx,n=n)
my.proj.migr <- PopProjWithMigration(Age=my.DF$Age,AgeGroup=my.DF$AgeGroup,bFx=my.DF$bFx,
                          sFx=my.DF$sFx,Ix=Ix,N0=my.DF$NFx,n=n)


## visualizing the difference between the two
## long data type
my.proj.nomigr.long <- my.proj.nomigr %>% 
  pivot_longer(-c(Age,AgeGroup),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
         year= (period-1)*5 + 1993,
         yearF = factor(year),
         type="no migration")
my.proj.migr.long <- my.proj.migr %>% 
  pivot_longer(-c(Age,AgeGroup),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
         year= (period-1)*5 + 1993,
         yearF = factor(year),
         type="with migration")
my.proj.long <- my.proj.nomigr.long %>% 
  bind_rows(my.proj.migr.long)

my.proj.long %>% 
  ggplot(aes(x=AgeGroup,y=population,fill=type)) +
  geom_bar(data = subset(my.proj.long, year %in% c(1998)),
           stat="identity",position = "dodge",color="black") +
  coord_flip() +
  ggtitle("Swedish female population") +
  theme_bw()


## ---- incorporating future assumptions

fx <- my.DF$Fx
plot(x,fx,t="l",lwd=2)
tfr.start <- 5*sum(fx)

tfr.seq <- seq(tfr.start,0.75,length.out=n)
plot(1:n,tfr.seq,t="l",lwd=2)

## create a matrix of time-specific fertility rates
FX <- matrix(NA,nrow = m,ncol=n)
FX[,1] <- fx

test <- tfr.seq[2]*fx/tfr.start
5*sum(test)

for (i in 2:n){
  FX[,i] <- tfr.seq[i]*fx/tfr.start
}

matplot(x,FX,lty=1,t="l",col=rainbow(n))

## make sure that the sum of FX is equal to our tfr series
5*apply(FX,2,sum)
tfr.seq


## make new pop proj function with time-specific fertiliy rates
PopProj.Fertility <- function(Age,AgeGroup,FX,sFx,N0,n,
                              srbF,LF0,l0){
  ## dimension of problem
  m <- length(N0)
  ## adjust the sFx and bFx vectors
  sFx <- sFx[-m]
  ## create the Leslie matrix
  L <- matrix(0,nrow=m,ncol=m)
  diag(L[-1,]) <- sFx
  L[m,m] <- sFx[m-1]
  ## create my output matrix
  N <- matrix(0,m,n+1)
  N[,1] <- N0
  i <- 1
  for (i in 1:n){
    bFx <- srbF*LF0/(2*l0) *(FX[-m,i] + sFx*FX[-1,i])
    L[1,1:length(bFx)] <- bFx
    N[,i+1] <- L%*%N[,i]
  }
  out <- cbind(data.frame(Age=Age,AgeGroup=AgeGroup),N)
  return(out)
}


n <- 20
my.proj.constant <- PopProj(Age=my.DF$Age,AgeGroup=my.DF$AgeGroup,bFx=my.DF$bFx,
                          sFx=my.DF$sFx,N0=my.DF$NFx,n=n)
my.proj.timespecific <- PopProj.Fertility(Age=my.DF$Age,
                            AgeGroup=my.DF$AgeGroup,FX=FX,
                            sFx=my.DF$sFx,N0=my.DF$NFx,n=n,
                            srbF=srbF,LF0=LF0,l0=l0)


## visualizing the difference between the two
## long data type
my.proj.constant.long <- my.proj.constant %>% 
  pivot_longer(-c(Age,AgeGroup),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
         year= (period-1)*5 + 1993,
         yearF = factor(year),
         type="constant")
my.proj.timespecific.long <- my.proj.timespecific %>% 
  pivot_longer(-c(Age,AgeGroup),names_to = "period",values_to = "population") %>% 
  mutate(period=as.numeric(period),
         year= (period-1)*5 + 1993,
         yearF = factor(year),
         type="time-varying")
my.proj.long <- my.proj.constant.long %>% 
  bind_rows(my.proj.timespecific.long)

my.proj.long %>% 
  ggplot(aes(x=AgeGroup,y=population,fill=type)) +
  geom_bar(data = subset(my.proj.long, year %in% c(2093)),
           stat="identity",position = "dodge",color="black") +
  coord_flip() +
  ggtitle("Swedish female population") +
  theme_bw()

## SHINY APP for dynamic visualization of your results
library(shiny)
ui <- fluidPage(
  sliderInput(inputId = "year", label = "Year", step = 5,
              value = min(my.proj.long$year), min = min(my.proj.long$year),
              max = max(my.proj.long$year)),
  column(12, plotOutput("plot_pyr1"))
)

server <- function(input, output){
  output$plot_pyr1 <- renderPlot({
    ## plotting pyramid
    ggplot(my.proj.long,aes(x=AgeGroup,y=population,fill=type)) +
      geom_bar(data = subset(my.proj.long, year == input$year),
               stat = "identity",position = "dodge",color = "black") +
      coord_flip() +
      theme_bw() +
      ggtitle(paste("Swedish female population, year",subset(my.proj.long, year == input$year)$year)) +
      scale_fill_manual(name = "Projection", values=c("#E69F00", "#56B4E9","#1C7C54")) +
      scale_y_continuous(limits = c(0, 350000), breaks = seq(0, 350000, 100000))
  })
}

shinyApp(ui = ui, server = server)


## alternative shiny app
library(plotly)
my.proj.long %>% 
  plot_ly(
    y=.~AgeGroup,
    x=.~population,
    frame=.~yearF,
    color=.~type,
    type="bar"
  )


## END

