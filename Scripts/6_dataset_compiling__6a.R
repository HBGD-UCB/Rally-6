

rm(list=ls())
library(caret)
library(tidyverse)

source("U:/R scripts/Wast_incidence_functions.R")
source("U:/Rally 6 India/R scripts/dataset_compiling_functions_6a.R")

setwd("U:/data/WastIncDatasets")


#------------------------------------------------
# Longitudinal dataset
#------------------------------------------------


compile_hbgdki_data(cum_inc=F, recoveryoutcome=F,
                    long.data=T,
                    data_location = "U:/data/",
                    file_location="U:/data/Compiled Datasets",
                    filename="6A_long.Rdata", 
                    rds=T,
                    suffix=NULL)



rm(list=ls())

#heatmap data
load("U:/data/Compiled Datasets/6A_long.Rdata")

d<-d[,-1]



df<-d %>% mutate(agemonth=factor(floor(AGEDAYS/30.25))) %>% 
  group_by(STUDYID, COUNTRY, agemonth) %>%
  summarize(
    N=n(), 
    meanWHZ=mean(WHZ, na.rm=T),
    meanHAZ=mean(HAZ, na.rm=T),
    wast=mean(as.numeric(WHZ <(-2))*100, na.rm=T),
    stunt=mean(as.numeric(HAZ <(-2))*100, na.rm=T),
    sevwast=mean(as.numeric(WHZ <(-3))*100, na.rm=T),
    sevstunt=mean(as.numeric(HAZ <(-3))*100, na.rm=T)
  )


save(df, file="U:/results/Rally 6A/heatmap_df.Rdata")













