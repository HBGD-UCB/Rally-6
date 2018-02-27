

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


d<-readRDS("zmrt.rds")
zmrt_inc<-WastIncCalc(d, washout=30)
zmrt_inc_table <- WastIncTable(zmrt_inc)
save(zmrt_inc, zmrt_inc_table, file="WastIncDatasets/zmrt_inc_30d.Rdata")




d<-readRDS("cmc.rds")
cmc_inc<-WastIncCalc(d, washout=30)
cmc_inc_table <- WastIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="WastIncDatasets/cmc_inc_30d.Rdata")



d<-readRDS("irc.rds")
irc_inc<-WastIncCalc(d, washout=30)
irc_inc_table <- WastIncTable(irc_inc)
save(irc_inc, irc_inc_table, file="WastIncDatasets/irc_inc_30d.Rdata")


d<-readRDS("tdc.rds")
tdc_inc<-WastIncCalc(d, washout=30)
tdc_inc_table <- WastIncTable(tdc_inc)
save(tdc_inc, tdc_inc_table, file="WastIncDatasets/tdc_inc_30d.Rdata")




d<-readRDS("mled.rds")
length(unique(d$COUNTRY))
mled_inc_bangladesh<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
mled_inc_brazil<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
mled_inc_india<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[3],])
mled_inc_nepal<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
mled_inc_peru<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
mled_inc_pakistan<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[6],])
mled_inc_southafrica<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[7],])
mled_inc_tanzania<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[8],])
mled_inc_table_bangladesh <- WastIncTable(mled_inc_bangladesh)
mled_inc_table_brazil <- WastIncTable(mled_inc_brazil)
mled_inc_table_india <- WastIncTable(mled_inc_india)
mled_inc_table_nepal <- WastIncTable(mled_inc_nepal)
mled_inc_table_peru <- WastIncTable(mled_inc_peru)
mled_inc_table_pakistan <- WastIncTable(mled_inc_pakistan)
mled_inc_table_southafrica <- WastIncTable(mled_inc_southafrica)
mled_inc_table_tanzania <- WastIncTable(mled_inc_tanzania)

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
     file="WastIncDatasets/mled_inc_30d.Rdata")
gc()



d<-readRDS("prvd.rds")
table(d$COUNTRY)
prvd_inc<-WastIncCalc(d, washout=30)
prvd_inc_table <- WastIncTable(prvd_inc)
save(prvd_inc, prvd_inc_table, file="WastIncDatasets/prvd_inc_30d.Rdata")

d<-readRDS("cmpf.rds")
table(d$COUNTRY)
cmpf_inc<-WastIncCalc(d, washout=30)
cmpf_inc_table <- WastIncTable(cmpf_inc)
save(cmpf_inc, cmpf_inc_table, file="WastIncDatasets/cmpf_inc_30d.Rdata")

d<-readRDS("fspp.rds")
table(d$COUNTRY)
fspp_inc<-WastIncCalc(d, washout=30)
fspp_inc_table <- WastIncTable(fspp_inc)
save(fspp_inc, fspp_inc_table, file="WastIncDatasets/fspp_inc_30d.Rdata")

d<-readRDS("cmc.rds")
table(d$COUNTRY)
cmc_inc<-WastIncCalc(d, washout=30)
cmc_inc_table <- WastIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="WastIncDatasets/cmc_inc_30d.Rdata")


d<-readRDS("eu.rds")
table(d$COUNTRY)
eu_inc<-WastIncCalc(d, washout=30)
eu_inc_table <- WastIncTable(eu_inc)
save(eu_inc, eu_inc_table, file="WastIncDatasets/eu_inc_30d.Rdata")



d<-readRDS("cort.rds")
unique(d$COUNTRY)
cort_inc_brazil<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
cort_inc_guatemala<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
cort_inc_india<-WastIncCalc(d, washout=30)
cort_inc_philippines<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
cort_inc_southafrica <- WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
cort_inc_table_brazil <- WastIncTable(cort_inc_brazil)
cort_inc_table_guatemala <- WastIncTable(cort_inc_guatemala)
cort_inc_table_india <- WastIncTable(cort_inc_india)
cort_inc_table_philippines <- WastIncTable(cort_inc_philippines)
cort_inc_table_southafrica <- WastIncTable(cort_inc_southafrica)
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
     file="WastIncDatasets/cort_inc_30d.Rdata")
gc()





d<-readRDS("dvds.rds")
table(d$COUNTRY)
dvds_inc<-WastIncCalc(d, washout=30)
dvds_inc_table <- WastIncTable(dvds_inc)
save(dvds_inc, dvds_inc_table, file="WastIncDatasets/dvds_inc_30d.Rdata")


d<-readRDS("vita.rds")
table(d$COUNTRY)
vita_inc<-WastIncCalc(d, washout=30)
vita_inc_table <- WastIncTable(vita_inc)
save(vita_inc, vita_inc_table, file="WastIncDatasets/vita_inc_30d.Rdata")

d<-readRDS("nvta.rds")
table(d$COUNTRY)
nvta_inc<-WastIncCalc(d, washout=30)
nvta_inc_table <- WastIncTable(nvta_inc)
save(nvta_inc, nvta_inc_table, file="WastIncDatasets/nvta_inc_30d.Rdata")

d<-readRDS("vb12.rds")
table(d$COUNTRY)
vb12_inc<-WastIncCalc(d, washout=30)
vb12_inc_table <- WastIncTable(vb12_inc)
save(vb12_inc, vb12_inc_table, file="WastIncDatasets/vb12_inc_30d.Rdata")

d<-readRDS("zlbw.rds")
table(d$COUNTRY)
zlbw_inc<-WastIncCalc(d, washout=30)
zlbw_inc_table <- WastIncTable(zlbw_inc)
save(zlbw_inc, zlbw_inc_table, file="WastIncDatasets/zlbw_inc_30d.Rdata")

d<-readRDS("ttbw.rds")
table(d$COUNTRY)
ttbw_inc<-WastIncCalc(d, washout=30)
ttbw_inc_table <- WastIncTable(ttbw_inc)
save(ttbw_inc, ttbw_inc_table, file="WastIncDatasets/ttbw_inc_30d.Rdata")

d<-readRDS("zsga.rds")
table(d$COUNTRY)
zsga_inc<-WastIncCalc(d, washout=30)
zsga_inc_table <- WastIncTable(zsga_inc)
save(zsga_inc, zsga_inc_table, file="WastIncDatasets/zsga_inc_30d.Rdata")

d<-readRDS("bts.rds")
table(d$COUNTRY)
bts_inc<-WastIncCalc(d, washout=30)
bts_inc_table <- WastIncTable(bts_inc)
save(bts_inc, bts_inc_table, file="WastIncDatasets/bts_inc_30d.Rdata")

d<-readRDS("nbsp.rds")
table(d$COUNTRY)
nbsp_inc<-WastIncCalc(d, washout=30)
nbsp_inc_table <- WastIncTable(nbsp_inc)
save(nbsp_inc, nbsp_inc_table, file="WastIncDatasets/nbsp_inc_30d.Rdata")

d<-readRDS("cbcy.rds")
table(d$COUNTRY)
cbcy_inc<-WastIncCalc(d, washout=30)
cbcy_inc_table <- WastIncTable(cbcy_inc)
save(cbcy_inc, cbcy_inc_table, file="WastIncDatasets/cbcy_inc_30d.Rdata")

d<-readRDS("zinf.rds")
table(d$COUNTRY)
zinf_inc<-WastIncCalc(d, washout=30)
zinf_inc_table <- WastIncTable(zinf_inc)
save(zinf_inc, zinf_inc_table, file="WastIncDatasets/zinf_inc_30d.Rdata")

d<-readRDS("imnc.rds")
table(d$COUNTRY)
imnc_inc<-WastIncCalc(d, washout=30)
imnc_inc_table <- WastIncTable(imnc_inc)
save(imnc_inc, imnc_inc_table, file="WastIncDatasets/imnc_inc_30d.Rdata")




#-------------------------------------------
# 90 day washout
#-------------------------------------------


d<-readRDS("zmrt.rds")
zmrt_inc<-WastIncCalc(d, washout=90)
zmrt_inc_table <- WastIncTable(zmrt_inc)
save(zmrt_inc, zmrt_inc_table, file="WastIncDatasets/zmrt_inc_90d.Rdata")




d<-readRDS("cmc.rds")
cmc_inc<-WastIncCalc(d, washout=90)
cmc_inc_table <- WastIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="WastIncDatasets/cmc_inc_90d.Rdata")



d<-readRDS("irc.rds")
irc_inc<-WastIncCalc(d, washout=90)
irc_inc_table <- WastIncTable(irc_inc)
save(irc_inc, irc_inc_table, file="WastIncDatasets/irc_inc_90d.Rdata")


d<-readRDS("tdc.rds")
tdc_inc<-WastIncCalc(d, washout=90)
tdc_inc_table <- WastIncTable(tdc_inc)
save(tdc_inc, tdc_inc_table, file="WastIncDatasets/tdc_inc_90d.Rdata")




d<-readRDS("mled.rds")
length(unique(d$COUNTRY))
mled_inc_bangladesh<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
mled_inc_brazil<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
mled_inc_india<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[3],])
mled_inc_nepal<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
mled_inc_peru<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
mled_inc_pakistan<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[6],])
mled_inc_southafrica<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[7],])
mled_inc_tanzania<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[8],])
mled_inc_table_bangladesh <- WastIncTable(mled_inc_bangladesh)
mled_inc_table_brazil <- WastIncTable(mled_inc_brazil)
mled_inc_table_india <- WastIncTable(mled_inc_india)
mled_inc_table_nepal <- WastIncTable(mled_inc_nepal)
mled_inc_table_peru <- WastIncTable(mled_inc_peru)
mled_inc_table_pakistan <- WastIncTable(mled_inc_pakistan)
mled_inc_table_southafrica <- WastIncTable(mled_inc_southafrica)
mled_inc_table_tanzania <- WastIncTable(mled_inc_tanzania)

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
     file="WastIncDatasets/mled_inc_90d.Rdata")
gc()



d<-readRDS("prvd.rds")
table(d$COUNTRY)
prvd_inc<-WastIncCalc(d, washout=90)
prvd_inc_table <- WastIncTable(prvd_inc)
save(prvd_inc, prvd_inc_table, file="WastIncDatasets/prvd_inc_90d.Rdata")

d<-readRDS("cmpf.rds")
table(d$COUNTRY)
cmpf_inc<-WastIncCalc(d, washout=90)
cmpf_inc_table <- WastIncTable(cmpf_inc)
save(cmpf_inc, cmpf_inc_table, file="WastIncDatasets/cmpf_inc_90d.Rdata")

d<-readRDS("fspp.rds")
table(d$COUNTRY)
fspp_inc<-WastIncCalc(d, washout=90)
fspp_inc_table <- WastIncTable(fspp_inc)
save(fspp_inc, fspp_inc_table, file="WastIncDatasets/fspp_inc_90d.Rdata")

d<-readRDS("cmc.rds")
table(d$COUNTRY)
cmc_inc<-WastIncCalc(d, washout=90)
cmc_inc_table <- WastIncTable(cmc_inc)
save(cmc_inc, cmc_inc_table, file="WastIncDatasets/cmc_inc_90d.Rdata")


d<-readRDS("eu.rds")
table(d$COUNTRY)
eu_inc<-WastIncCalc(d, washout=90)
eu_inc_table <- WastIncTable(eu_inc)
save(eu_inc, eu_inc_table, file="WastIncDatasets/eu_inc_90d.Rdata")



d<-readRDS("cort.rds")
unique(d$COUNTRY)
cort_inc_brazil<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[1],])
cort_inc_guatemala<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[2],])
cort_inc_india<-WastIncCalc(d, washout=90)
cort_inc_philippines<-WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[4],])
cort_inc_southafrica <- WastIncCalc(d[d$COUNTRY==unique(d$COUNTRY)[5],])
cort_inc_table_brazil <- WastIncTable(cort_inc_brazil)
cort_inc_table_guatemala <- WastIncTable(cort_inc_guatemala)
cort_inc_table_india <- WastIncTable(cort_inc_india)
cort_inc_table_philippines <- WastIncTable(cort_inc_philippines)
cort_inc_table_southafrica <- WastIncTable(cort_inc_southafrica)
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
     file="WastIncDatasets/cort_inc_90d.Rdata")
gc()





d<-readRDS("dvds.rds")
table(d$COUNTRY)
dvds_inc<-WastIncCalc(d, washout=90)
dvds_inc_table <- WastIncTable(dvds_inc)
save(dvds_inc, dvds_inc_table, file="WastIncDatasets/dvds_inc_90d.Rdata")


d<-readRDS("vita.rds")
table(d$COUNTRY)
vita_inc<-WastIncCalc(d, washout=90)
vita_inc_table <- WastIncTable(vita_inc)
save(vita_inc, vita_inc_table, file="WastIncDatasets/vita_inc_90d.Rdata")

d<-readRDS("nvta.rds")
table(d$COUNTRY)
nvta_inc<-WastIncCalc(d, washout=90)
nvta_inc_table <- WastIncTable(nvta_inc)
save(nvta_inc, nvta_inc_table, file="WastIncDatasets/nvta_inc_90d.Rdata")

d<-readRDS("vb12.rds")
table(d$COUNTRY)
vb12_inc<-WastIncCalc(d, washout=90)
vb12_inc_table <- WastIncTable(vb12_inc)
save(vb12_inc, vb12_inc_table, file="WastIncDatasets/vb12_inc_90d.Rdata")

d<-readRDS("zlbw.rds")
table(d$COUNTRY)
zlbw_inc<-WastIncCalc(d, washout=90)
zlbw_inc_table <- WastIncTable(zlbw_inc)
save(zlbw_inc, zlbw_inc_table, file="WastIncDatasets/zlbw_inc_90d.Rdata")

d<-readRDS("ttbw.rds")
table(d$COUNTRY)
ttbw_inc<-WastIncCalc(d, washout=90)
ttbw_inc_table <- WastIncTable(ttbw_inc)
save(ttbw_inc, ttbw_inc_table, file="WastIncDatasets/ttbw_inc_90d.Rdata")

d<-readRDS("zsga.rds")
table(d$COUNTRY)
zsga_inc<-WastIncCalc(d, washout=90)
zsga_inc_table <- WastIncTable(zsga_inc)
save(zsga_inc, zsga_inc_table, file="WastIncDatasets/zsga_inc_90d.Rdata")

d<-readRDS("bts.rds")
table(d$COUNTRY)
bts_inc<-WastIncCalc(d, washout=90)
bts_inc_table <- WastIncTable(bts_inc)
save(bts_inc, bts_inc_table, file="WastIncDatasets/bts_inc_90d.Rdata")

d<-readRDS("nbsp.rds")
table(d$COUNTRY)
nbsp_inc<-WastIncCalc(d, washout=90)
nbsp_inc_table <- WastIncTable(nbsp_inc)
save(nbsp_inc, nbsp_inc_table, file="WastIncDatasets/nbsp_inc_90d.Rdata")

d<-readRDS("cbcy.rds")
table(d$COUNTRY)
cbcy_inc<-WastIncCalc(d, washout=90)
cbcy_inc_table <- WastIncTable(cbcy_inc)
save(cbcy_inc, cbcy_inc_table, file="WastIncDatasets/cbcy_inc_90d.Rdata")

d<-readRDS("zinf.rds")
table(d$COUNTRY)
zinf_inc<-WastIncCalc(d, washout=90)
zinf_inc_table <- WastIncTable(zinf_inc)
save(zinf_inc, zinf_inc_table, file="WastIncDatasets/zinf_inc_90d.Rdata")

d<-readRDS("imnc.rds")
table(d$COUNTRY)
imnc_inc<-WastIncCalc(d, washout=90)
imnc_inc_table <- WastIncTable(imnc_inc)
save(imnc_inc, imnc_inc_table, file="WastIncDatasets/imnc_inc_90d.Rdata")
