

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
source("Stunt_incidence_functions.R")

setwd("U:/data")


d<-readRDS("zmrt.rds")
zmrt_inc<-StuntIncCalc(d, washout=60)
zmrt_inc_table <- StuntIncTable(zmrt_inc)
save(zmrt_inc, zmrt_inc_table, file="StuntIncDatasets/zmrt_st_inc.Rdata")




d<-readRDS("cmc.rds")
cmc_inc<-StuntIncCalc(d, washout=60)
cmc_inc_table <- StuntIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="StuntIncDatasets/cmc_st_inc.Rdata")



d<-readRDS("irc.rds")
irc_inc<-StuntIncCalc(d, washout=60)
irc_inc_table <- StuntIncTable(irc_inc)
save(irc_inc, irc_inc_table, file="StuntIncDatasets/irc_st_inc.Rdata")


d<-readRDS("tdc.rds")
tdc_inc<-StuntIncCalc(d, washout=60)
tdc_inc_table <- StuntIncTable(tdc_inc)
save(tdc_inc, tdc_inc_table, file="StuntIncDatasets/tdc_st_inc.Rdata")




d<-readRDS("mled.rds")
length(unique(d$COUNTRY))
mled_inc_bangladesh<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
mled_inc_brazil<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
mled_inc_india<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[3],])
mled_inc_nepal<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
mled_inc_peru<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
mled_inc_pakistan<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[6],])
mled_inc_southafrica<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[7],])
mled_inc_tanzania<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[8],])
mled_inc_table_bangladesh <- StuntIncTable(mled_inc_bangladesh)
mled_inc_table_brazil <- StuntIncTable(mled_inc_brazil)
mled_inc_table_india <- StuntIncTable(mled_inc_india)
mled_inc_table_nepal <- StuntIncTable(mled_inc_nepal)
mled_inc_table_peru <- StuntIncTable(mled_inc_peru)
mled_inc_table_pakistan <- StuntIncTable(mled_inc_pakistan)
mled_inc_table_southafrica <- StuntIncTable(mled_inc_southafrica)
mled_inc_table_tanzania <- StuntIncTable(mled_inc_tanzania)

save(mled_inc_bangladesh,
     mled_inc_brazil,
     mled_inc_india,
     mled_inc_nepal,
     mled_inc_peru,
     mled_inc_pakistan,
     mled_inc_southafrica,
     mled_inc_tanzania,
     mled_inc_table_bangladesh,
     mled_inc_table_brazil,
     mled_inc_table_india,
     mled_inc_table_nepal,
     mled_inc_table_peru,
     mled_inc_table_pakistan,
     mled_inc_table_southafrica,
     mled_inc_table_tanzania,
     file="StuntIncDatasets/mled_st_inc.Rdata")
gc()



d<-readRDS("prvd.rds")
table(d$COUNTRY)
prvd_inc<-StuntIncCalc(d, washout=60)
prvd_inc_table <- StuntIncTable(prvd_inc)
save(prvd_inc, prvd_inc_table, file="StuntIncDatasets/prvd_st_inc.Rdata")

d<-readRDS("cmpf.rds")
table(d$COUNTRY)
cmpf_inc<-StuntIncCalc(d, washout=60)
cmpf_inc_table <- StuntIncTable(cmpf_inc)
save(cmpf_inc, cmpf_inc_table, file="StuntIncDatasets/cmpf_st_inc.Rdata")

d<-readRDS("fspp.rds")
table(d$COUNTRY)
fspp_inc<-StuntIncCalc(d, washout=60)
fspp_inc_table <- StuntIncTable(fspp_inc)
save(fspp_inc, fspp_inc_table, file="StuntIncDatasets/fspp_st_inc.Rdata")

d<-readRDS("cmc.rds")
table(d$COUNTRY)
cmc_inc<-StuntIncCalc(d, washout=60)
cmc_inc_table <- StuntIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="StuntIncDatasets/cmc_st_inc.Rdata")


d<-readRDS("eu.rds")
table(d$COUNTRY)
eu_inc<-StuntIncCalc(d, washout=60)
eu_inc_table <- StuntIncTable(eu_inc)
save(eu_inc, eu_inc_table, file="StuntIncDatasets/eu_st_inc.Rdata")



d<-readRDS("cort.rds")
unique(d$COUNTRY)
cort_inc_brazil<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
cort_inc_guatemala<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
cort_inc_india<-StuntIncCalc(d, washout=60)
cort_inc_philippines<-StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
cort_inc_southafrica <- StuntIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
cort_inc_table_brazil <- StuntIncTable(cort_inc_brazil)
cort_inc_table_guatemala <- StuntIncTable(cort_inc_guatemala)
cort_inc_table_india <- StuntIncTable(cort_inc_india)
cort_inc_table_philippines <- StuntIncTable(cort_inc_philippines)
cort_inc_table_southafrica <- StuntIncTable(cort_inc_southafrica)
save(cort_inc_brazil,
     cort_inc_guatemala,
     cort_inc_india,
     cort_inc_philippines,
     cort_inc_southafrica,
     cort_inc_table_brazil,
     cort_inc_table_guatemala,
     cort_inc_table_india,
     cort_inc_table_philippines,
     cort_inc_table_southafrica,
     file="StuntIncDatasets/cort_st_inc.Rdata")
gc()





d<-readRDS("dvds.rds")
table(d$COUNTRY)
dvds_inc<-StuntIncCalc(d, washout=60)
dvds_inc_table <- StuntIncTable(dvds_inc)
save(dvds_inc, dvds_inc_table, file="StuntIncDatasets/dvds_st_inc.Rdata")


d<-readRDS("vita.rds")
table(d$COUNTRY)
vita_inc<-StuntIncCalc(d, washout=60)
vita_inc_table <- StuntIncTable(vita_inc)
save(vita_inc, vita_inc_table, file="StuntIncDatasets/vita_st_inc.Rdata")

d<-readRDS("nvta.rds")
table(d$COUNTRY)
nvta_inc<-StuntIncCalc(d, washout=60)
nvta_inc_table <- StuntIncTable(nvta_inc)
save(nvta_inc, nvta_inc_table, file="StuntIncDatasets/nvta_st_inc.Rdata")

d<-readRDS("vb12.rds")
table(d$COUNTRY)
vb12_inc<-StuntIncCalc(d, washout=60)
vb12_inc_table <- StuntIncTable(vb12_inc)
save(vb12_inc, vb12_inc_table, file="StuntIncDatasets/vb12_st_inc.Rdata")

d<-readRDS("zlbw.rds")
table(d$COUNTRY)
zlbw_inc<-StuntIncCalc(d, washout=60)
zlbw_inc_table <- StuntIncTable(zlbw_inc)
save(zlbw_inc, zlbw_inc_table, file="StuntIncDatasets/zlbw_st_inc.Rdata")

d<-readRDS("ttbw.rds")
table(d$COUNTRY)
ttbw_inc<-StuntIncCalc(d, washout=60)
ttbw_inc_table <- StuntIncTable(ttbw_inc)
save(ttbw_inc, ttbw_inc_table, file="StuntIncDatasets/ttbw_st_inc.Rdata")

d<-readRDS("zsga.rds")
table(d$COUNTRY)
zsga_inc<-StuntIncCalc(d, washout=60)
zsga_inc_table <- StuntIncTable(zsga_inc)
save(zsga_inc, zsga_inc_table, file="StuntIncDatasets/zsga_st_inc.Rdata")

d<-readRDS("bts.rds")
table(d$COUNTRY)
bts_inc<-StuntIncCalc(d, washout=60)
bts_inc_table <- StuntIncTable(bts_inc)
save(bts_inc, bts_inc_table, file="StuntIncDatasets/bts_st_inc.Rdata")

d<-readRDS("nbsp.rds")
table(d$COUNTRY)
nbsp_inc<-StuntIncCalc(d, washout=60)
nbsp_inc_table <- StuntIncTable(nbsp_inc)
save(nbsp_inc, nbsp_inc_table, file="StuntIncDatasets/nbsp_st_inc.Rdata")

d<-readRDS("cbcy.rds")
table(d$COUNTRY)
cbcy_inc<-StuntIncCalc(d, washout=60)
cbcy_inc_table <- StuntIncTable(cbcy_inc)
save(cbcy_inc, cbcy_inc_table, file="StuntIncDatasets/cbcy_st_inc.Rdata")

d<-readRDS("zinf.rds")
table(d$COUNTRY)
zinf_inc<-StuntIncCalc(d, washout=60)
zinf_inc_table <- StuntIncTable(zinf_inc)
save(zinf_inc, zinf_inc_table, file="StuntIncDatasets/zinf_st_inc.Rdata")

d<-readRDS("imnc.rds")
table(d$COUNTRY)
imnc_inc<-StuntIncCalc(d, washout=60)
imnc_inc_table <- StuntIncTable(imnc_inc)
save(imnc_inc, imnc_inc_table, file="StuntIncDatasets/imnc_st_inc.Rdata")



