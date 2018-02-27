

rm(list=ls())
library(dplyr)
library(tidyverse)
library(caret)
library(MASS)
library(reshape2)
library(zoo)
library(epitools)
library(binom)
theme_set(theme_bw())

setwd("U:/R scripts")
source("Wast_incidence_functions.R")

setwd("U:/data")

d<-readRDS("eu.rds")
table(d$COUNTRY)
eu_inc<-WastIncCalc(d, dropBornWasted=T)
eu_inc_table <- WastIncTable(eu_inc)
save(eu_inc, eu_inc_table, file="WastIncDatasets/eu_inc_NoBirthInc.Rdata")


d<-readRDS("dvds.rds")
table(d$COUNTRY)
dvds_inc<-WastIncCalc(d, dropBornWasted=T)
dvds_inc_table <- WastIncTable(dvds_inc)
save(dvds_inc, dvds_inc_table, file="WastIncDatasets/dvds_inc_NoBirthInc.Rdata")


d<-readRDS("vita.rds")
table(d$COUNTRY)
vita_inc<-WastIncCalc(d, dropBornWasted=T)
vita_inc_table <- WastIncTable(vita_inc)
save(vita_inc, vita_inc_table, file="WastIncDatasets/vita_inc_NoBirthInc.Rdata")

d<-readRDS("nvta.rds")
table(d$COUNTRY)
nvta_inc<-WastIncCalc(d, dropBornWasted=T)
nvta_inc_table <- WastIncTable(nvta_inc)
save(nvta_inc, nvta_inc_table, file="WastIncDatasets/nvta_inc_NoBirthInc.Rdata")

d<-readRDS("vb12.rds")
table(d$COUNTRY)
vb12_inc<-WastIncCalc(d, dropBornWasted=T)
vb12_inc_table <- WastIncTable(vb12_inc)
save(vb12_inc, vb12_inc_table, file="WastIncDatasets/vb12_inc_NoBirthInc.Rdata")

d<-readRDS("zlbw.rds")
table(d$COUNTRY)
zlbw_inc<-WastIncCalc(d, dropBornWasted=T)
zlbw_inc_table <- WastIncTable(zlbw_inc)
save(zlbw_inc, zlbw_inc_table, file="WastIncDatasets/zlbw_inc_NoBirthInc.Rdata")

d<-readRDS("ttbw.rds")
table(d$COUNTRY)
ttbw_inc<-WastIncCalc(d, dropBornWasted=T)
ttbw_inc_table <- WastIncTable(ttbw_inc)
save(ttbw_inc, ttbw_inc_table, file="WastIncDatasets/ttbw_inc_NoBirthInc.Rdata")

d<-readRDS("zsga.rds")
table(d$COUNTRY)
zsga_inc<-WastIncCalc(d, dropBornWasted=T)
zsga_inc_table <- WastIncTable(zsga_inc)
save(zsga_inc, zsga_inc_table, file="WastIncDatasets/zsga_inc_NoBirthInc.Rdata")

d<-readRDS("bts.rds")
table(d$COUNTRY)
bts_inc<-WastIncCalc(d, dropBornWasted=T)
bts_inc_table <- WastIncTable(bts_inc)
save(bts_inc, bts_inc_table, file="WastIncDatasets/bts_inc_NoBirthInc.Rdata")

d<-readRDS("nbsp.rds")
table(d$COUNTRY)
nbsp_inc<-WastIncCalc(d, dropBornWasted=T)
nbsp_inc_table <- WastIncTable(nbsp_inc)
save(nbsp_inc, nbsp_inc_table, file="WastIncDatasets/nbsp_inc_NoBirthInc.Rdata")

d<-readRDS("cbcy.rds")
table(d$COUNTRY)
cbcy_inc<-WastIncCalc(d, dropBornWasted=T)
cbcy_inc_table <- WastIncTable(cbcy_inc)
save(cbcy_inc, cbcy_inc_table, file="WastIncDatasets/cbcy_inc_NoBirthInc.Rdata")

d<-readRDS("zinf.rds")
table(d$COUNTRY)
zinf_inc<-WastIncCalc(d, dropBornWasted=T)
zinf_inc_table <- WastIncTable(zinf_inc)
save(zinf_inc, zinf_inc_table, file="WastIncDatasets/zinf_inc_NoBirthInc.Rdata")

d<-readRDS("imnc.rds")
table(d$COUNTRY)
imnc_inc<-WastIncCalc(d, dropBornWasted=T)
imnc_inc_table <- WastIncTable(imnc_inc)
save(imnc_inc, imnc_inc_table, file="WastIncDatasets/imnc_inc_NoBirthInc.Rdata")




