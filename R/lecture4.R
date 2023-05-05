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

