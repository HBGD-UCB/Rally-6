#-----------------------------------
# Stunting analysis
# Objective 1a
# Import data, subset to relevant variables
#-----------------------------------
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(data.table)

setwd("U:/data/GHAP_data/")
load("U:/data/Stunting/Full-compiled-data/compiled_HAZ_dataset.RData")

#--------------------------------------------
# Subset to relevant variables
#--------------------------------------------
colnames(d)=tolower(colnames(d))
d <- d %>% select(studyid, subjid, country, tr, agedays, haz)
d <- d%>% filter(country=="INDIA")

nrow(d)

#--------------------------------------------
# drop acutely ill enrollment studies
#--------------------------------------------
d <- d %>% filter(studyid!="ki1000301-DIVIDS")

#--------------------------------------------
# Add in TDC
#--------------------------------------------

tdc <- readRDS("U:/data/tdc.rds")
colnames(tdc) <- tolower(colnames(tdc))
tdc <- tdc %>% subset(., select=c(studyid, subjid, country, agedays, haz))
tdc$tr <- NA
class(tdc$subjid) <- "integer64"

d <- bind_rows(d, tdc)

#--------------------------------------------
# drop unrealistic HAZ
#--------------------------------------------
nrow(d)
d = filter(d,haz >= -6 & haz <=6)
nrow(d)

#--------------------------------------------
# order data, create measurement id
#--------------------------------------------
d <- d %>% 
  arrange(studyid,subjid,agedays) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid,agedays) %>%
  # create id for measurement within person
  mutate(measid=seq_along(subjid)) 

# count number of studies
length(names(table(d$studyid)))

# table of studies
table(d$studyid)
table(d$studyid,d$country)




save(d,file="U:/Data/Stunting/stunting_data_6A.RData")



#tabulate children and observations
df <- d%>% filter(country=="INDIA") %>% filter(haz > -6 & haz < 6) %>% filter(agedays < 25*30.4167)

length(unique(d$studyid))
length(unique(paste0(d$studyid,"-",d$subjid)))
dim(d)

