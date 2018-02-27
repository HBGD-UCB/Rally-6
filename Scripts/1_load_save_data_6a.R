









rm(list=ls())
library(dplyr)
library(tidyr)
library(caret)
library(ghap)

setwd("U:/data")
set_git_base_path("U:/git")
get_git_base_path()

d <- get_study_list_anthro()
meta<-d[grepl("IND",d$country),]
save(meta, file="U:/Rally 6 India/Results/metadata6a.Rdata")

#Create function to fill in static covariates
impute_static_vars<-function(d){
  
  
  varlist<-c("STUDYID","SUBJID", "shortid",
             "SITEID","SEXN",
             "SEX","AGEDAYS",
             "CTRYCD","COUNTRY",
             "WTKG","WAZ",
             "HAZ","BAZ",
             "BMI","WHZ",
             "LENCM","LATITUDE",
             "LONGITUD","HTCM",
             "SUBJIDO","CITYTOWN",
             "ELEVATN","BIRTHWT",
             "MUACCM","MAGE",
             "MUAZ","MEDUCYRS",
             "HCIRCM","GAGEBRTH",
             "ARMCD","ARM",
             "SANITATN","NPERSON",
             "MHTCM",
             "BIRTHLEN","FEEDINGN",
             "FEEDING","TV",
             "BRTHYR","REGCTRY",
             "RADIO",
             "SUMEP","SUMDIAR",
             "SUMDAYS","PCTDIAR",
             "DLVLOC",
             "BRTHWEEK","FLOOR",
             "DELIVERY",
             "AGEDTH","FEDUCYRS",
             "PARITY","NROOMS",
             "CAR","COOKFUEL",
             "WALL",
             "VISITNUM","VISIT",
             "BICYCLE","MWORK",
             "H2OSRCP","ELEC",
             "FRIG","GAGEDAYS",
             "DEAD","NCHLDLT5",
             "ROOF",
             "MMARITN","MMARIT",
             "H2OTRTP","MWTKG",
             "MCYCLE","INCTOT",
             "APGAR5","FWORK",
             "H2OSRC","GOAT",
             "BIRTHHC","DURBRST",
             "MBMI","SEWING",
             "NCHLD","FAGE",
             "NLCHILD","PHONE",
             "COW",
             "APGAR1",
             "SMOKED","AGLAND",
             "COUGHFL","CHCIRCM",
             "LLPHONE","BED",
             "CLTHCAB","COOKPLAC",
             "NADULT","FHTCM",
             "CHICKEN","VISITDY",
             "SESN","SES",
             "SHEEP","MULTBRTH",
             "WATCH","CHAIR",
             "TABLE","STRATUMN",
             "STRATUM","EDUCCRGV",
             "RACE",
             "APGAR10","M_WTKG",
             "BFEDFL","DOGS",
             "CATTLE","H2OSRCC",
             "MRACE","METHNIC",
             "OWNHOME",
             "FAN",
             "ANIMALS","CART",
             "WASHMAC","SMOKAMT",
             "NFCHILD","NMCHILD",
             "H2OSRCB",
             "ETHNIC","M_HTCM",
             "FWTKG","FHOUSEH",
             "CATS",
             "AGLANDSZ",
             "FOODDFCT",
             "NUMHCVIS",
             "PREGOUT","RICKSHAW",
             "MATTRESS","H2OAVAIL",
             "ANTPTNUM",
             "ANTPT","MAGE1ST",
             "NSLEEP","KITCHDSC",
             "M_BMI",
             "PREECLMP",
             "SOAP","EDUCHH",
             "MSMKSTAT",
             "RODENTS",
             "VITD",
             "H2OSRCS",
             "FOWL","STOVE",
             "MHOUSEH","H2OFREQ",
             "H2OUNTRT","WASHCCHD",
             "WASHCOOK","WASHNURS",
             "R_MUAZ",
             "BEDNET","H2OSRCK",
             "FRACE",
             "SINK","SOFA",
             "NBEDROOM",
             "FSMKSTAT",
             "COHORTN","COHORT",
             "EDUCFAM",
             "INCMOM","PIGS",
             "CASTE",
             "MBRTHWT",
             "MTRIBE","FTRIBE",
             "WTKGM","HTCMM",
             "BMIM","WAZM",
             "HAZM","WHZM",
             "BAZM",
             "SMOKSTAT","SMOKYRS",
             "DIARFL",
             #WBK variables
             "STUDYID", "SUBJID",  "CLUSTID", "SITEID",  "HHID",    "SEXN",    "SEX",     "ARMCD",   "ARM",     "STRATUMN","STRATUM",
             "COHORTN", "COHORT",  "BRTHYR",  "BRTHWEEK","MHTCM",   "MEDUCYRS","CTRYCD",  "COUNTRY", "BICYCLE", "CAR",     "CHICKEN",
             "COW",     "DOGS",    "ELEC",    "GOAT",    "MOBILE",  "MCYCLE",  "HUNGER",  "RADIO",   "STOVE",   "TV",      "H2OTIME",
             "NCOMP",   "WATCH",   "FLOOR",   "COOKFUEL","ROOF",    "SANITATN","WALL",    "CHLDLT18","AGEDAYS", "VISITNUM","VISIT",  
             "WTKG",    "LENCM",   "BMI",     "HCIRCM",  "WAZ",     "HAZ",     "WHZ",     "BAZ",     "HCAZ",    "FREECHL", "IMPRLAT",
             "WATSOAP", "LNSN",    "LNSP",    "FCSAFDIS",
             #additional diarrhea vars
             "DIARFL","SUMDIAR",  
             "DIARDAYS", "CSUMDIAR",
             "DIAREPS" , "DIARBFEN" ,
             "DIARRHOEA","DIARRHOEA_NEONATAL")
  
  
  dynamicvars<-c(
    "AGEDAYS",
    "WTKG","WAZ",
    "HAZ","BAZ",
    "BMI","WHZ",
    "LENCM","HTCM",
    "HCIRCM",
    "SUMEP","SUMDIAR",
    "SUMDAYS","PCTDIAR",
    "VISITNUM","VISIT",
    "DEAD","DURBRST",
    "COUGHFL","CHCIRCM",
    "VISITDY",
    "ANTPTNUM",
    "FREECHL", "IMPRLAT",
    "WATSOAP", "LNSN",    "LNSP",    "FCSAFDIS",
    #additional diarrhea vars
    "DIARFL","SUMDIAR",  
    "DIARDAYS", "CSUMDIAR",
    "DIAREPS" , "DIARBFEN" ,
    "DIARRHOEA","DIARRHOEA_NEONATAL")
  
  
  
  
  
  
  
  staticvars<-varlist[-which(varlist %in% dynamicvars)]
  
  study.d <- d %>% filter(!is.na(AGEDAYS) & !is.na(HAZ))
  
  if(!is.null(varlist)){
    study.d<-study.d[,which(colnames(study.d) %in% varlist)]
  }
  study.d <- apply(study.d, 2, as.character)
  study.d<-as.data.frame(study.d)
  
  
  #Set "" and other missing codes to missing
  for(i in 1:ncol(study.d)){
    study.d[,i]<-ifelse(study.d[,i]=="",NA,as.character(study.d[,i]))
  } 
  
  #seperate out dynamic variables
  study.d.varying<-study.d[,which(colnames(study.d) %in% c("SUBJID", "AGEDAYS", dynamicvars))]
  
  #fill in missing static variables
  study.d.static<-study.d[,which(colnames(study.d) %in% c("SUBJID", "AGEDAYS", staticvars))]
  
  study.d.static<-study.d.static %>%  
    group_by(SUBJID) %>%
    do(fill(.,everything())) %>% 
    do(fill(.,everything(), .direction = 'up')) 
  
  study.d <- merge(study.d.static, study.d.varying, by=c("SUBJID","AGEDAYS"))
  
  study.d$AGEDAYS<-as.numeric(as.character(study.d$AGEDAYS))
  study.d$HAZ<-as.numeric(as.character(study.d$HAZ))
  study.d$WHZ<-as.numeric(as.character(study.d$WHZ))
  study.d$SUBJID<-as.numeric(as.character(study.d$SUBJID))
  
  return(study.d)
}





setwd("U:/data/")


d<-use_study("VITAMIN-A") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="vita.rds")
rm(d)


d<-use_study("NEOVITA") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="nvta.rds")
rm(d)


d<-use_study("Vitamin-B12") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="vb12.rds")
rm(d)


d<-use_study("LBW") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="zlbw.rds")
rm(d)

d<-use_study("TTBW") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ttbw_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="ttbw.rds")
rm(d)

d<-use_study("ZincSGA") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="zsga.rds")
rm(d)



d<-use_study("BtS") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="bts.rds")
rm(d)


d<-use_study("NBSP") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="nbsp.rds")
rm(d)


d<-use_study("CordBloodCytometry") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="ncry_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="cbcy.rds")
rm(d)


d<-use_study("zinf") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="zinf_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="zinf.rds")
rm(d)


d<-use_study("dvds") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="dvds_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="dvds.rds")
rm(d)


d<-use_study("znmort")        
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="zmrt_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="zmrt.rds")
rm(d) 



d<-use_study("imnci") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="imnc_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="imnc.rds")
rm(d)

d<-use_study("cmc_v_bcs_2002") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="cmc_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="cmc.rds")
rm(d)  

d<-use_study("cohorts")   
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="cort_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="cort.rds")
rm(d) 


d<-use_study("eu") 
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="eu_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="eu.rds")
rm(d)   


d<-use_study("irc")            
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="irc_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="irc.rds")
rm(d)  

#d<-use_study("mal_ed") #Mal ED not loading
#Temp fix
d<-read.csv("U:/data/MALED-201501/adam/trellisdat.csv")
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="mled_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="mled.rds")
rm(d) 


d<-use_study("sas_compfeed")   
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="cmpf_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="cmpf.rds")
rm(d) 

d<-use_study("sas_foodsuppl")  
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="fspp_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="fspp.rds")
rm(d) 

d<-use_study("tdc")            
colnames(d)<- toupper(colnames(d))
saveRDS(d, file="tdc_raw.rds")
d <- impute_static_vars(d)
saveRDS(d, file="tdc.rds")
rm(d)  






