## -------------------------------------------- ##
## EDSD 2022-2023: Population projections
## Lecture 2
## The cohort component method
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
load("data/dta.swe.1993.Rdata")

## start by inspecting our data
head(dta.swe)

## example of lead function
DF.example <- as_tibble(dta.swe) %>% 
  select(Age,NFx,LFx) %>% 
  mutate(LFx5=lead(LFx,n=3),
         sFx=LFx5/LFx,
         NFxlag=lag(NFx),
         dummy=ifelse(test = Age==10,
                      yes =  1,
                      no  =  0),
         dummy2=case_when(
           Age==10 ~ 1,
           Age==20 ~ 1,
           Age==30 ~ 3,
           TRUE   ~ 0
         ))

## compute sFx
my.DF <- as_tibble(dta.swe) %>% 
  mutate(sFx=lead(LFx)/LFx,
         NFx5 = lag(NFx*sFx),
         sFx = ifelse(test = Age == 80,
                     yes  = lead(LFx)/(LFx+lead(LFx)),
                     no   = sFx),
         NFx5 = ifelse(test = Age == 85,
                       yes  = (NFx+lag(NFx))*lag(sFx),
                       no   = NFx5))

## checking
my.DF$sFx[1]
my.DF$LFx[2]/my.DF$LFx[1]
my.DF$NFx5[2]
my.DF$NFx[1]*my.DF$sFx[1]


## adjusting the first age group
SRB <- 1.05
srbF <- 1/(1+SRB)
LF0 <- dta.swe$LFx[1]
l0 <- 1e5

my.DF <- my.DF %>% 
  mutate(bFx=srbF*LF0/(2*l0) *(Fx + sFx*lead(Fx)),
         Bx = Fx*5*(NFx+NFx5)/2,
         NFx5 = ifelse(test = Age == 0,
                       yes  = srbF*LF0/(5*l0) * sum(Bx,na.rm = T),
                       no   = NFx5))

my.DF$NFx5[1]

## make sure that the formulas are the same
sum(my.DF$NFx * my.DF$bFx,na.rm = T)

## population pyramid

## long data type
my.DF.long <- my.DF %>% 
  select(AgeGroup,NFx,NFx5) %>% 
  rename('1993'=NFx,'1998'=NFx5) %>% 
  pivot_longer(-AgeGroup,names_to = "year",values_to = "population")

## pyramid
my.DF.long %>% 
  ggplot(aes(x=AgeGroup,y=population,fill=year)) +
  geom_bar(stat="identity",position = "dodge",color="black") +
  coord_flip() +
  ggtitle("Swedish female population")
  
  




