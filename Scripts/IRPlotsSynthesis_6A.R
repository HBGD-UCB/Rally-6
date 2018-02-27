
# load packages
rm(list=ls())
library(tidyverse)

setwd("U:/data/WastIncDatasets")


#-------------------------------------
# Counting birth incidence
#-------------------------------------
load("cmc_inc.Rdata")
load("eu_inc.Rdata")
load("irc_inc.Rdata")
load("tdc_inc.Rdata")
load("cort_inc.Rdata")
load("mled_inc.Rdata")
load("dvds_inc.Rdata")
load("vita_inc.Rdata")
load("vb12_inc.Rdata")
load("zlbw_inc.Rdata")
load("zsga_inc.Rdata")
load("zinf_inc.Rdata")
load("cmpf_inc.Rdata")
load("fspp_inc.Rdata")
load("zmrt_inc.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  dvds_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)

#Make list of summary tables
tablelist <- list(
  cmc_inc_table, 
  eu_inc_table, 
  irc_inc_table,
  tdc_inc_table,
  cort_inc_table_india,            
  mled_inc_table_india,
  dvds_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table,
  cmpf_inc_table,
  fspp_inc_table,
  zmrt_inc_table)

#Define the corresponding study names
studynames <- c(
  "CMC Vellore Birth Cohort 2002",
  "Zn Supp for Diarrheal Pneuomonia", 
  "IRC (Vellore Crypto Study)",
  "TDC (Vellore Bottled Water BC)",
  "COHORTS Study: India",
  "MAL-ED Study: India",
  "DIVIDS-1 with 5 year Follow-up",
  "EPI Linked Vita A",
  "Vitamin B12 supplementation",
  "Zn Supp Trial in LBW",
  "Zinc supplementation in infants born small for gestational age",
  "Zn Treatment for Diarrhea",
  "SAS-CompFeed: Optimal Infant Feeding Practices",
  "SAS-FoodSuppl: Randomized Food Supplementation Trial Ages 4-12 Months",
  "Zn Supp and Infant Mortality"
)

d<-NULL
for(i in 1:length(datasets)){
  cat(datasets[[i]]$STUDYID[1], " ", datasets[[i]]$COUNTRY[1])
  temp<-cbind(rep(datasets[[i]]$STUDYID[1], nrow(tablelist[[i]]$means)), rep(datasets[[i]]$COUNTRY[1],nrow(tablelist[[i]]$means)), tablelist[[i]]$means)
  d<-rbind(d,temp)
}  
colnames(d)[1:2]<-c("STUDYID","COUNTRY")


#Save summary statistics list to create pooled estimates

setwd("U:/Rally 6 India/Results")
save(d, file="descriptive_epi_mean_monthly_cohorts_6A.Rdata")



#-------------------------------------
# No birth incidence
#-------------------------------------

setwd("U:/data/WastIncDatasets")

load("cmc_inc_NoBirthInc.Rdata")
load("eu_inc_NoBirthInc.Rdata")
load("irc_inc_NoBirthInc.Rdata")
load("tdc_inc_NoBirthInc.Rdata")
load("cort_inc_NoBirthInc.Rdata")
load("mled_inc_NoBirthInc.Rdata")
load("dvds_inc_NoBirthInc.Rdata")
load("vita_inc_NoBirthInc.Rdata")
load("vb12_inc_NoBirthInc.Rdata")
load("zlbw_inc_NoBirthInc.Rdata")
load("zsga_inc_NoBirthInc.Rdata")
load("zinf_inc_NoBirthInc.Rdata")
load("cmpf_inc_NoBirthInc.Rdata")
load("fspp_inc_NoBirthInc.Rdata")
load("zmrt_inc_NoBirthInc.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  dvds_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)

#Make list of summary tables
tablelist <- list(
  cmc_inc_table, 
  eu_inc_table, 
  irc_inc_table,
  tdc_inc_table,
  cort_inc_table_india,            
  mled_inc_table_india,
  dvds_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table,
  cmpf_inc_table,
  fspp_inc_table,
  zmrt_inc_table)

#Define the corresponding study names
studynames <- c(
  "CMC Vellore Birth Cohort 2002",
  "Zn Supp for Diarrheal Pneuomonia", 
  "IRC (Vellore Crypto Study)",
  "TDC (Vellore Bottled Water BC)",
  "COHORTS Study: India",
  "MAL-ED Study: India",
  "DIVIDS-1 with 5 year Follow-up",
  "EPI Linked Vita A",
  "Vitamin B12 supplementation",
  "Zn Supp Trial in LBW",
  "Zinc supplementation in infants born small for gestational age",
  "Zn Treatment for Diarrhea",
  "SAS-CompFeed: Optimal Infant Feeding Practices",
  "SAS-FoodSuppl: Randomized Food Supplementation Trial Ages 4-12 Months",
  "Zn Supp and Infant Mortality"
)

d<-NULL
for(i in 1:length(datasets)){
  cat(datasets[[i]]$STUDYID[1], " ", datasets[[i]]$COUNTRY[1])
  temp<-cbind(rep(datasets[[i]]$STUDYID[1], nrow(tablelist[[i]]$means)), rep(datasets[[i]]$COUNTRY[1],nrow(tablelist[[i]]$means)), tablelist[[i]]$means)
  d<-rbind(d,temp)
}  
colnames(d)[1:2]<-c("STUDYID","COUNTRY")


setwd("U:/Rally 6 India/Results")
save(d, file="descriptive_epi_mean_monthly_cohorts_noBW_6A.Rdata")



#-------------------------------------
# 30 day washout
#-------------------------------------
setwd("U:/data/WastIncDatasets")



load("cmc_inc_30d.Rdata")
load("eu_inc_30d.Rdata")
load("irc_inc_30d.Rdata")
load("tdc_inc_30d.Rdata")
load("cort_inc_30d.Rdata")
load("mled_inc_30d.Rdata")
load("dvds_inc_30d.Rdata")
load("vita_inc_30d.Rdata")
load("vb12_inc_30d.Rdata")
load("zlbw_inc_30d.Rdata")
load("zsga_inc_30d.Rdata")
load("zinf_inc_30d.Rdata")
load("cmpf_inc_30d.Rdata")
load("fspp_inc_30d.Rdata")
load("zmrt_inc_30d.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  dvds_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)

#Make list of summary tables
tablelist <- list(
  cmc_inc_table, 
  eu_inc_table, 
  irc_inc_table,
  tdc_inc_table,
  cort_inc_table_india,            
  mled_inc_table_india,
  dvds_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table,
  cmpf_inc_table,
  fspp_inc_table,
  zmrt_inc_table)

#Define the corresponding study names
studynames <- c(
  "CMC Vellore Birth Cohort 2002",
  "Zn Supp for Diarrheal Pneuomonia", 
  "IRC (Vellore Crypto Study)",
  "TDC (Vellore Bottled Water BC)",
  "COHORTS Study: India",
  "MAL-ED Study: India",
  "DIVIDS-1 with 5 year Follow-up",
  "EPI Linked Vita A",
  "Vitamin B12 supplementation",
  "Zn Supp Trial in LBW",
  "Zinc supplementation in infants born small for gestational age",
  "Zn Treatment for Diarrhea",
  "SAS-CompFeed: Optimal Infant Feeding Practices",
  "SAS-FoodSuppl: Randomized Food Supplementation Trial Ages 4-12 Months",
  "Zn Supp and Infant Mortality"
)

d<-NULL
for(i in 1:length(datasets)){
  cat(datasets[[i]]$STUDYID[1], " ", datasets[[i]]$COUNTRY[1])
  temp<-cbind(rep(datasets[[i]]$STUDYID[1], nrow(tablelist[[i]]$means)), rep(datasets[[i]]$COUNTRY[1],nrow(tablelist[[i]]$means)), tablelist[[i]]$means)
  d<-rbind(d,temp)
}  
colnames(d)[1:2]<-c("STUDYID","COUNTRY")

setwd("U:/Rally 6 India/Results")
save(d, file="descriptive_epi_mean_monthly_cohorts_30d_6A.Rdata")


#-------------------------------------
# 90 day washout
#-------------------------------------

setwd("U:/data/WastIncDatasets")


load("cmc_inc_90d.Rdata")
load("eu_inc_90d.Rdata")
load("irc_inc_90d.Rdata")
load("tdc_inc_90d.Rdata")
load("cort_inc_90d.Rdata")
load("mled_inc_90d.Rdata")
load("dvds_inc_90d.Rdata")
load("vita_inc_90d.Rdata")
load("vb12_inc_90d.Rdata")
load("zlbw_inc_90d.Rdata")
load("zsga_inc_90d.Rdata")
load("zinf_inc_90d.Rdata")
load("cmpf_inc_90d.Rdata")
load("fspp_inc_90d.Rdata")
load("zmrt_inc_90d.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  dvds_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)

#Make list of summary tables
tablelist <- list(
  cmc_inc_table, 
  eu_inc_table, 
  irc_inc_table,
  tdc_inc_table,
  cort_inc_table_india,            
  mled_inc_table_india,
  dvds_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table,
  cmpf_inc_table,
  fspp_inc_table,
  zmrt_inc_table)

#Define the corresponding study names
studynames <- c(
  "CMC Vellore Birth Cohort 2002",
  "Zn Supp for Diarrheal Pneuomonia", 
  "IRC (Vellore Crypto Study)",
  "TDC (Vellore Bottled Water BC)",
  "COHORTS Study: India",
  "MAL-ED Study: India",
  "DIVIDS-1 with 5 year Follow-up",
  "EPI Linked Vita A",
  "Vitamin B12 supplementation",
  "Zn Supp Trial in LBW",
  "Zinc supplementation in infants born small for gestational age",
  "Zn Treatment for Diarrhea",
  "SAS-CompFeed: Optimal Infant Feeding Practices",
  "SAS-FoodSuppl: Randomized Food Supplementation Trial Ages 4-12 Months",
  "Zn Supp and Infant Mortality"
)

d<-NULL
for(i in 1:length(datasets)){
  cat(datasets[[i]]$STUDYID[1], " ", datasets[[i]]$COUNTRY[1])
  temp<-cbind(rep(datasets[[i]]$STUDYID[1], nrow(tablelist[[i]]$means)), rep(datasets[[i]]$COUNTRY[1],nrow(tablelist[[i]]$means)), tablelist[[i]]$means)
  d<-rbind(d,temp)
}  
colnames(d)[1:2]<-c("STUDYID","COUNTRY")

setwd("U:/Rally 6 India/Results")
save(d, file="descriptive_epi_mean_monthly_cohorts_90d_6A.Rdata")


#-------------------------------------
# Stunting incidence
#-------------------------------------

setwd("U:/data/StuntIncDatasets")


load("cmc_st_inc.Rdata")
load("eu_st_inc.Rdata")
load("irc_st_inc.Rdata")
load("tdc_st_inc.Rdata")
load("cort_st_inc.Rdata")
load("mled_st_inc.Rdata")
load("dvds_st_inc.Rdata")
load("vita_st_inc.Rdata")
load("vb12_st_inc.Rdata")
load("zlbw_st_inc.Rdata")
load("zsga_st_inc.Rdata")
load("zinf_st_inc.Rdata")
load("cmpf_st_inc.Rdata")
load("fspp_st_inc.Rdata")
load("zmrt_st_inc.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  dvds_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)

#Make list of summary tables
tablelist <- list(
  cmc_inc_table, 
  eu_inc_table, 
  irc_inc_table,
  tdc_inc_table,
  cort_inc_table_india,            
  mled_inc_table_india,
  dvds_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table,
  cmpf_inc_table,
  fspp_inc_table,
  zmrt_inc_table)

#Define the corresponding study names
studynames <- c(
  "CMC Vellore Birth Cohort 2002",
  "Zn Supp for Diarrheal Pneuomonia", 
  "IRC (Vellore Crypto Study)",
  "TDC (Vellore Bottled Water BC)",
  "COHORTS Study: India",
  "MAL-ED Study: India",
  "DIVIDS-1 with 5 year Follow-up",
  "EPI Linked Vita A",
  "Vitamin B12 supplementation",
  "Zn Supp Trial in LBW",
  "Zinc supplementation in infants born small for gestational age",
  "Zn Treatment for Diarrhea",
  "SAS-CompFeed: Optimal Infant Feeding Practices",
  "SAS-FoodSuppl: Randomized Food Supplementation Trial Ages 4-12 Months",
  "Zn Supp and Infant Mortality"
)

d<-NULL
for(i in 1:length(datasets)){
  cat(datasets[[i]]$STUDYID[1], " ", datasets[[i]]$COUNTRY[1])
  temp<-cbind(rep(datasets[[i]]$STUDYID[1], nrow(tablelist[[i]]$means)), rep(datasets[[i]]$COUNTRY[1],nrow(tablelist[[i]]$means)), tablelist[[i]]$means)
  d<-rbind(d,temp)
}  
colnames(d)[1:2]<-c("STUDYID","COUNTRY")





setwd("U:/Rally 6 India/Results")
save(d, file="descriptive_epi_mean_monthly_cohorts_stunting_6A.Rdata")


