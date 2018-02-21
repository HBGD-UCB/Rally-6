





#http://www.reed.edu/data-at-reed/software/R/markdown_multiple_reports.html
#https://stackoverflow.com/questions/30422008/r-knitr-pdf-is-there-a-posssibility-to-automatically-save-pdf-reports-generate

# References for automation 
# http://www.r-bloggers.com/how-to-source-an-r-script-automatically-on-a-mac-using-automator-and-ical/
# http://www.engadget.com/2013/03/18/triggering-applescripts-from-calendar-alerts-in-mountain-lion/

# File 1: Should be an R-Script 
# contains a loop that iteratively calls an Rmarkdown file (i.e. File 2)

# load packages
rm(list=ls())
library(knitr)
library(markdown)
library(rmarkdown)
library(ggplot2)

setwd("U:/R scripts/")
source("Markdown_plot_functions.R")

setwd("U:/data/WastIncDatasets")


load("cort_inc.Rdata")
load("eu_inc.Rdata")
load("irc_inc.Rdata")
load("mled_inc.Rdata")
load("cmpf_inc.Rdata")
load("fspp_inc.Rdata")
load("tdc_inc.Rdata")
load("dvds_inc.Rdata")
load("cmc_inc.Rdata")
load("tdc_inc.Rdata")
load("zmrt_inc.Rdata")
load("vita_inc.Rdata")
load("vb12_inc.Rdata")
load("zlbw_inc.Rdata")
load("ttbw_inc.Rdata")
load("zsga_inc.Rdata")
load("bts_inc.Rdata")
load("cbcy_inc.Rdata")
load("zinf_inc.Rdata")
load("imnc_inc.Rdata")



#Define set of short ID's I care about
shortids <- c(
  "CORT-India", 
  "EU", "IRC",
  "MLED-India",
  "CMPF","FSPP","TDC","DVDS",
  "CMC", "ZMRT", "VITA",  "VB12", "ZLBW", "ZSGA",
  "ZINF"
)



#List of datasets
datasets <- list(
  cort_inc_india,
  eu_inc,
  irc_inc,
  mled_inc_india,
  cmpf_inc,
  fspp_inc,
  tdc_inc,
  dvds_inc,
  cmc_inc,
  zmrt_inc,
  vita_inc,
  vb12_inc,
  zlbw_inc,
  zsga_inc,
  zinf_inc)




#Make list of summary tables
tablelist<- list(
  cort_inc_table_india,
  eu_inc_table,
  irc_inc_table,
  mled_inc_table_india,
  cmpf_inc_table,
  fspp_inc_table,
  tdc_inc_table,
  dvds_inc_table,
  cmc_inc_table,
  zmrt_inc_table,
  vita_inc_table,
  vb12_inc_table,
  zlbw_inc_table,
  zsga_inc_table,
  zinf_inc_table)








#Define the corresponding study names
studynames <- c(
  "cort COHORTS India",
  "eu Zn Supp for Diarrheal Pneuomonia",
  "irc Vellore Crypto Study",
  "mled MAL-ED Study India",
  "cmpf Optimal Infant Feeding Practices",
  "fspp Randomized Food Supplementation Trial Ages 4-12 Months",
  "tdc Transmission Dynamics of Cryptosporidial Infections",
  "dvds DIVIDS-1 with 5 year Follow-up",
  "cmc CMC Vellore Birth Cohort 2002" ,
  "zmrt Zn Supp and Infant Mortality",
  "vita EPI Linked Vita A",
  "vb12 Vitamin B12 supplementation",
  "zlbw Zn Supp Trial in LBW",
  "zsga Zinc supplementation in infants born small for gestational age",
  "zinf Zn Treatment for Diarrhea"
)




# render allows the rmarkdown to still have access to the global environment, 
#so I just need to define the variables I want read in within the for-loop
for(i in 1:length(datasets)){
  study <- studynames[i]
  d <- datasets[[i]]
  tables <- tablelist[[i]]
  
  
  
  try(
    rmarkdown::render('U:/Rally 6 India/R scripts/MarkdownTemplate_6A.Rmd',  
                      output_file =  paste("report_",shortids[i],'_', Sys.Date(), ".html", sep=''), 
                      output_dir = 'U:/Rally 6 India/Results/Reports')
  )}



