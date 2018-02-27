

#rm(list=ls())
#library(dplyr)
#library(tidyr)
#library(SuperLearner)
#library(caret)


#setwd("U:/results/Stunting")
#load("st_GHAPstudies_raw.Rdata")


#-------------------------------------------------------------------
# Impute missingness function
#-------------------------------------------------------------------


Mode<-function(x){
  ux = unique(x)
  tab = tabulate(match(x, ux))
  ux[tab == max(tab)]
}

missingness_indicators<-function (data, prefix = "miss_", remove_constant = T, remove_collinear = T, 
                                  skip_vars = c(), verbose = F) {
  indicators = sapply(data[, !colnames(data) %in% skip_vars], 
                      FUN = function(col) as.numeric(is.na(col)))
  colnames(indicators) = paste0(prefix, colnames(data)[!colnames(data) %in% 
                                                         skip_vars])
  if (remove_constant) {
    col_means = colMeans(indicators)
    if (verbose) {
      num_removed = sum(col_means %in% c(0, 1))
      if (num_removed > 0) {
        cat("Removing", num_removed, "indicators that are constant.\n")
      }
    }
    indicators = indicators[, !col_means %in% c(0, 1), drop = F]
  }
  if (remove_collinear) {
    linear_combos = caret::findLinearCombos(indicators)
    remove_columns = linear_combos$remove
    if (length(linear_combos$remove) > 0) {
      if (verbose) {
        cat("Removing", length(linear_combos$remove), 
            "indicators due to collinearity:\n")
        cat(paste0(colnames(indicators)[linear_combos$remove], 
                   collapse = ", "), "\n")
      }
      indicators = indicators[, -linear_combos$remove, 
                              drop = F]
    }
  }
  return(indicators)
}

impute_missing_values<-function(data, type = "standard", add_indicators = T, prefix = "miss_", skip_vars = c(), verbose = F){
  missing_indicators = NULL
  new_data = data
  results = list(type = type, add_indicators = add_indicators, 
                 skip_vars = skip_vars, prefix = prefix)
  if (type == "standard") {
    if (verbose) {
      cat("Running standard imputation.\n")
    }
    preprocess = NA
    impute_values = vector("list", ncol(data))
    names(impute_values) = colnames(data)
    for (i in 1:ncol(data)) {
      nas = sum(is.na(data[[i]]))
      col_class = class(data[[i]])
      if (col_class == "factor") {
        impute_value = Mode(data[[i]])[1]
      }
      else if (col_class %in% c("integer", "numeric")) {
        impute_value = median(data[[i]], na.rm = T)
      }
      else {
        warning(paste(colnames(data)[i], "should be numeric or factor type. But its class is", 
                      col_class))
      }
      impute_values[[i]] = impute_value
      if (nas == 0 || names(data)[i] %in% skip_vars) {
        next
      }
      else if (nas == nrow(data)) {
        if (verbose) {
          cat("Note: skipping", colnames(data)[i], "because all values are NA.\n")
        }
        next
      }
      else {
        new_data[is.na(data[[i]]), i] = impute_value
      }
    }
    results$impute_values = impute_values
  }
  else if (type == "knn") {
    impute_info = caret::preProcess(new_data, method = c("knnImpute"))
    new_data = predict(impute_info, new_data)
    results$impute_info = impute_info
  }
  if (add_indicators) {
    missing_indicators = missingness_indicators(data, prefix = prefix, 
                                                verbose = verbose)
    if (verbose) {
      cat("Indicators added:", ncol(missing_indicators), 
          "\n")
    }
    new_data = cbind(new_data, missing_indicators)
  }
  results$data = new_data
  return(results)
}


impute_missing_hbgdki <- function(d,  outcomes = c("stunt", "sevstunt", "HTCM", "LENCM",
                                                   "WHZ", "HAZ",  "wast","sevwast","period_length","sevwast_inc",
                                                   "wast_inc","sevwast_rec","wast_rec","wast_risk",
                                                   "sevwast_risk","wast_rec_risk","sevwast_rec_risk","wasting_episode",
                                                   "born_wast_inc","episode_ID","incident_age","maxage",
                                                   "duration","wasting_duration","born_sevwast_inc","sevwasting_episode",
                                                   "sev_episode_ID","sevwasting_duration", "wast_rec90d","wast_rec60d",
                                                   "wast_rec30d","sevwast_inc90d","sevwast_inc60d","sevwast_inc30d")){
  #require(ck37r)
  
  outcome.df <- d[,which(colnames(d) %in% c("STUDYID", "SUBJID", outcomes))]
  d <- d[,which(!(colnames(d) %in% outcomes))]
  
  #Impute overall median to create missingness indicator columns
  missdf <- impute_missing_values(as.data.frame(d), type = "standard", add_indicators = T,prefix = "miss_", verbose = F)$data
  
  #add in missing indicator columns
  dim(d)
  dim(missdf)
  d <- data.frame(d, missdf[,(ncol(d)+1):ncol(missdf)])
  dim(d)
  
  
  #Add in missing columns with no variance for risk factors with no missingness
  d$miss_SEX<-0
  d$miss_region<-0
  d$miss_month<-0
  d$miss_enrolstunt<-0
  
  #Fill in study specific median for studies where a variable was measured
  d<- d %>% group_by(STUDYID, COUNTRY) %>% 
    do(as.data.frame(impute_missing_values(as.data.frame(.), type = "standard", add_indicators = F, verbose = F)$data))
  
  df <- merge(outcome.df, d, by=c("STUDYID","SUBJID"))
  
  return(df)
}




#-------------------------------------------------------------------
# Slice observation at age of interest and bind rows with other ghap studies
#-------------------------------------------------------------------


bindGHAP<-function(study, varlist=NULL, dynamicvars=NULL, d, 
                   age=24*30.25, agerange=c(12*30.25,36*30.25), minage=0, maxage=24*30.25,
                   cum_inc=F, recoveryoutcome=F, long.data=F, rds=T, noBW=F, suffix="_inc", filesuffix=NULL, country=NULL){
  
  if(rds==T){
    study.d<-readRDS(paste0(study,".rds"))
  }else{
    object=paste0(study,suffix)
    if(is.null(filesuffix)){
      load(paste0(object,".Rdata"))
    }else{
      load(paste0(study,filesuffix,".Rdata"))
    }
    if(!is.null(country)){
      object=paste0(study,suffix,"_",country)
      
    }
    assign("study.d",get(object))
    rm(object)
  }

  
  if("WHZ" %in% colnames(study.d)){
    
    study.d <- apply(study.d, 2, as.character)
    study.d <- as.data.frame(study.d)
    study.d$AGEDAYS<-as.numeric(as.character(study.d$AGEDAYS))
    study.d$HAZ<-as.numeric(as.character(study.d$HAZ))
    study.d$WHZ<-as.numeric(as.character(study.d$WHZ))
    study.d$SUBJID<-as.numeric(as.character(study.d$SUBJID))
    
    
    study.d$age <- age
    study.d$agerange1 <- agerange[1]
    study.d$agerange2 <- agerange[2]
    study.d$minage <- minage
    study.d$maxage <- maxage
    
    
    #Find anthro at enrollment
    enrol_anthro <-study.d %>% group_by(SUBJID) %>% 
      arrange(AGEDAYS) %>% 
      subset(., select=c(COUNTRY,SUBJID,HAZ, WHZ)) %>%
      slice(1) %>% #Find enrollment observation 
      rename(enrolHAZ=HAZ, enrolWHZ=WHZ) %>%
      ungroup
    
    if(recoveryoutcome==T){
      cat(study,": recovery outcome\n")
      
      
      study.d <-study.d %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
        filter(AGEDAYS > minage & AGEDAYS <= maxage) %>%
        filter(wast_inc==1) %>%
        slice(1) %>% 
        ungroup
    }else{
      
      if(long.data==T){
        cat(study,": long data\n")
        
        study.d <- study.d[,which(colnames(study.d) %in% varlist)]
        
        study.d <-study.d %>% group_by(SUBJID) %>% 
          arrange(AGEDAYS) %>%
          mutate(anywast=ifelse(sum(WHZ < (-2), na.rm=T)>0,1,0), anysevwast=ifelse(sum(WHZ < (-2), na.rm=T)>0,1,0),
                 anystunt=ifelse(sum(HAZ < (-2), na.rm=T)>0,1,0), anysevstunt=ifelse(sum(HAZ < (-2), na.rm=T)>0,1,0)) %>%
          ungroup
      }else{
        
        if(cum_inc==T){
          if(noBW==T){
            cat(study,": CI no birth wasting\n")
            
            temp <-study.d %>% group_by(SUBJID) %>% 
              filter(AGEDAYS > minage & AGEDAYS <= maxage) %>%
              arrange(AGEDAYS) %>%
              mutate(wast_inc=ifelse(sum(wast_inc, na.rm=T)>0,1,0), sevwast_inc=ifelse(sum(sevwast_inc, na.rm=T)>0,1,0)) %>%
              slice(1) %>% 
              ungroup
          }else{
            cat(study,": CI\n")
            
            study.d <-study.d %>% group_by(SUBJID) %>% 
              filter(AGEDAYS > minage & AGEDAYS <= maxage) %>%
              arrange(AGEDAYS) %>%
              mutate(wast_inc=ifelse(sum(WHZ < (-2), na.rm=T)>0,1,0), sevwast_inc=ifelse(sum(WHZ < (-3), na.rm=T)>0,1,0)) %>%
              slice(1) %>% 
              ungroup 
          }
        }else{
          cat(study,": prevalence\n")
          
          #find age closest to selected age
          study.d <-study.d %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
            slice(which.min(abs(AGEDAYS - age))) %>% #Find observation with age closest to 2 years
            filter(AGEDAYS > agerange1 & AGEDAYS < agerange2) %>% #Select ages between 1 and 3 years
            ungroup
        }
      }
    }
    #merge in enrollment anthro
    study.d <- left_join(study.d, enrol_anthro, by=c("COUNTRY","SUBJID"))
    
  }else{
    study.d <- apply(study.d, 2, as.character)
    study.d <- as.data.frame(study.d)
    study.d$AGEDAYS<-as.numeric(as.character(study.d$AGEDAYS))
    study.d$HAZ<-NA
    study.d$WHZ<-NA
    study.d$SUBJID<-as.numeric(as.character(study.d$SUBJID))
  }
  
  if(nrow(study.d)>0){
    if(long.data==F){
    study.d <- subset(study.d, select = -c(age, agerange1, agerange2, minage, maxage))
    }
    d<- bind_rows(d, study.d , .id = "studyid")
  }
  rm(study.d)
  return(d)
}


#-------------------------------------------------------------------
# Wrapper function to bind included studies together
#-------------------------------------------------------------------

compile_hbgdki_data <- function(age=24*30.25, agerange=c(12*30.25, 36*30.25), minage=0, maxage=24*30.25,
                                cum_inc=F, recoveryoutcome=F,
                                data_location="U:/data/WastIncDatasets",
                                file_location="U:/results",
                                filename, long.data=F,
                                rds=T, noBW=F,
                                suffix="_inc",
                                filesuffix=NULL){
  
  vars<-c("STUDYID","SUBJID", "shortid",
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
  
  
  dynvars<-c(
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
  
  
  setwd(data_location)
  
  
  d<-NULL 
  d<-bindGHAP(study="mled", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, country="india", long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="cort", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, country="india", long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="cmc", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW)  
  d<-bindGHAP(study="eu", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW)   
  d<-bindGHAP(study="irc", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW)  
  d<-bindGHAP(study="tdc", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW)  
  d<-bindGHAP(study="fspp", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="zmrt", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="cmpf", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="vita", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="nvta", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="vb12", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="zlbw", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="ttbw", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="zsga", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="bts", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="nbsp", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="cbcy", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="zinf", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="dvds", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  d<-bindGHAP(study="imnc", varlist=vars, dynamicvars=dynvars, d=d,age=age, agerange=agerange, minage=minage, maxage=maxage, cum_inc=cum_inc, recoveryoutcome=recoveryoutcome ,rds=rds, suffix=suffix, long.data=long.data, filesuffix=filesuffix, noBW=noBW) 
  
  
  # #clean covariates and impute missingness
  if(long.data==F){
    d <- clean_covariates_hbgdki(d)
    d <- impute_missing_hbgdki(d)
  }
  
  setwd(file_location)
  save(d, file=filename)
}



#-------------------------------------------------------------------
# Clean covariate function
#-------------------------------------------------------------------

clean_covariates_hbgdki <- function(d){
  
  d<-d[,-1]
  d$STUDYID<-as.factor(d$STUDYID)
  
  
  #Drop Mal-ED Pakistan
  d<- d[!(d$STUDYID=="MAL-ED" & d$COUNTRY=="PAKISTAN"),]
  
  
  
  #Drop unrealistic values of HAZ
  dim(d)
  d<-d %>% filter(HAZ>-5 & HAZ<5)
  dim(d)
  summary(d$HAZ)
  
  
  #Calculate wasting
  d$wast<-ifelse(d$WHZ< (-2), 1,0)
  d$wast[is.na(d$WHZ)]<-NA
  
  d$sevwast<-ifelse(d$WHZ< (-3), 1,0)
  d$sevwast[is.na(d$WHZ)]<-NA
  
  #Calculate stunting
  d$stunt<-ifelse(d$HAZ< (-2), 1,0)
  d$stunt[is.na(d$HAZ)]<-NA
  
  d$sevstunt<-ifelse(d$HAZ< (-3), 1,0)
  d$sevstunt[is.na(d$HAZ)]<-NA
  
  
  
  
  
  #################################################
  # Child characteristics
  #################################################
  
  #------------------------------------------------
  # 1) Sex
  #------------------------------------------------
  table(d$SEX)
  d$SEX[d$SEX=="male"]<-"Male"
  d$SEX[d$SEX=="female"]<-"Female"
  d$SEX<-as.factor(d$SEX)
  
  #------------------------------------------------
  # 2) Gestational age at birth
  #------------------------------------------------
  
  d$GAGEBRTH <- as.numeric(d$GAGEBRTH)
  summary(d$GAGEBRTH)
  
  #------------------------------------------------
  # 3) Birthweight
  #------------------------------------------------
  
  d$BIRTHWT <- as.numeric(d$BIRTHWT)
  d$BIRTHLEN <- as.numeric(d$BIRTHLEN)
  
  summary(d$BIRTHWT)
  summary(d$BIRTHLEN)
  
  #Drop unrealistic values
  table(d$BIRTHWT[d$BIRTHWT>3400])
  d$BIRTHWT[d$BIRTHWT>7000] <- NA
  
  table(d$BIRTHLEN[d$BIRTHLEN>51])
  d$BIRTHLEN[d$BIRTHLEN>70] <- NA
  
  #------------------------------------------------
  # 4) Stunted at enrollment or birth
  #------------------------------------------------
  
  
  d$enrolstunt<-ifelse(d$enrolHAZ<(-2),1,0)
  d$enrolstunt[is.na(d$enrolHAZ)]<-NA
  
  table(d$enrolstunt)
  
  
  
  #------------------------------------------------
  # 5) Birth order
  #------------------------------------------------
  
  d$PARITY <- as.numeric(d$PARITY)
  
  
  table(d$PARITY)
  table(d$STUDYID,d$PARITY)
  
  #harmonize so 1 is firstborn for all
  d <-d %>% group_by(STUDYID) %>% mutate(PARITY = if(any(min(PARITY, na.rm=T)==0))  PARITY+1 else  PARITY)
  
  
  d$birthorder <- NA
  d$birthorder[d$PARITY==1] <- "firstborn"
  d$birthorder[d$PARITY==2] <- "secondborn"
  d$birthorder[d$PARITY>2] <- "thirdborn+"
  table(d$birthorder)
  
  #Note: Where is WASH Benefits?
  
  #------------------------------------------------
  # 6) Calendar month of birth
  #------------------------------------------------
  
  #"BRTHWEEK"
  d$BRTHWEEK <- as.numeric(d$BRTHWEEK)
  
  table(d$BRTHWEEK)
  
  d$birthmonth <- round(d$BRTHWEEK/4.25)
  d$birthmonth[d$birthmonth==0]<-1
  table(d$birthmonth)
  
  
  #------------------------------------------------
  # 7) Delivery location (hospital)
  #------------------------------------------------
  
  
  table(d$DLVLOC)
  d$homedelivery <- ifelse(d$DLVLOC=="Home" | d$DLVLOC=="Maternity home",1,0)
  d$homedelivery[is.na(d$DLVLOC)] <- NA
  table(d$homedelivery)
  
  #------------------------------------------------
  # 8) Delivery method (C-section vs. vaginal)
  #------------------------------------------------
  
  d$DELIVERY <- as.factor(d$DELIVERY)
  
  table(d$DELIVERY)
  
  d$vagbirth <- ifelse(d$DELIVERY=="Vaginal" | 
                         d$DELIVERY=="Breech Vaginal" |
                         d$DELIVERY=="Forceps" |
                         d$DELIVERY=="Normal Vaginal" |
                         d$DELIVERY=="Normal vaginal birth" |
                         d$DELIVERY=="Normal Viginal" |
                         d$DELIVERY=="Spontaneous vaginal" |
                         d$DELIVERY=="Spontaneous Vaginal Delivery" |
                         d$DELIVERY=="Suction" |
                         d$DELIVERY=="Vacuum" |
                         d$DELIVERY=="Vacuum Extraction" ,1,0)
  d$vagbirth[is.na(d$DELIVERY)] <- NA
  table(d$vagbirth)
  
  
  #################################################
  # Maternal characteristics
  #################################################
  
  #------------------------------------------------
  # 9) Height
  #------------------------------------------------
  #Combine different variable names for maternal height
  d[is.na(d$MHTCM),"MHTCM"]<-d[is.na(d$MHTCM),"M_HTCM"]
  
  d$MHTCM <- as.numeric(d$MHTCM)
  
  
  #------------------------------------------------
  # 10) Weight
  #------------------------------------------------
  
  
  
  #Combine different variable names for maternal weight and BMI
  d[is.na(d$MWTKG),"MWTKG"]<-d[is.na(d$MWTKG),"M_WTKG"]
  d[is.na(d$MBMI),"MBMI"]<-d[is.na(d$MBMI),"M_BMI"]
  d$MWTKG <- as.numeric(d$MWTKG)
  d$MBMI <- as.numeric(d$MBMI)
  
  summary(d$MWTKG)
  summary(d$MBMI)
  
  
  #------------------------------------------------
  # 11) Age
  #------------------------------------------------
  
  d$MAGE <- as.numeric(d$MAGE)
  summary(d$MAGE)
  
  #Figure out which study does not have age classified in years
  d %>% group_by(STUDYID) %>% summarize(mn=mean(MAGE,na.rm=T)) %>% filter(!is.na(mn))
  
  d$MAGE[d$STUDYID=="ki1119695-PROBIT"] <- d$MAGE[d$STUDYID=="ki1119695-PROBIT"]/365.25
  
  d %>% group_by(STUDYID) %>% summarize(mn=mean(MAGE,na.rm=T)) %>% filter(!is.na(mn))
  
  
  #------------------------------------------------
  # 12) Education
  #------------------------------------------------
  d$MEDUCYRS<-as.numeric(as.character(d$MEDUCYRS))
  
  #------------------------------------------------
  # 13) Marital status
  #------------------------------------------------
  
  #FHOUSEH Father lives in household
  table(d$FHOUSEH)
  table(d$MMARIT)
  
  d$single <- ifelse(d$FHOUSEH=="No, permanently elsewhere",1,0)
  d$single[is.na(d$FHOUSEH)] <- NA
  
  d$single <- ifelse(d$MMARIT=="Common law" | d$MMARIT=="Married",0,1)
  d$single[is.na(d$MMARIT) & is.na(d$FHOUSEH)] <- NA
  table(d$single)
  
  #------------------------------------------------
  # 14) Breastfeeding practices 
  #------------------------------------------------
  d$breastfeeding<-"Exclusively breastfed"
  d$breastfeeding[d$FEEDING=="Mixture breast/formula fed" | d$FEEDING!=""]<-"Mixed"
  d$breastfeeding[d$FEEDING=="Exclusively formula fed" | d$FEEDING=="Bottle"]<-"Exclusively formula fed"
  d$breastfeeding[d$FEEDING=="Unknown"]<-"Missing"
  d$breastfeeding<-factor(d$breastfeeding, levels=c("Exclusively breastfed", "Mixed", "Exclusively formula fed"))
  table(d$breastfeeding)
  
  
  #################################################
  # Paternal characteristics
  #################################################
  
  #------------------------------------------------
  # 15) Height
  #------------------------------------------------
  
  #Combine different variable names for paternal height
  d$FHTCM <- as.numeric(d$FHTCM)
  summary(d$FHTCM)
  
  #------------------------------------------------
  # 16) Age
  #------------------------------------------------
  d$FAGE <- as.numeric(d$FAGE)
  summary(d$FAGE)
  
  #Figure out which study does not have age classified in years
  d %>% group_by(STUDYID) %>% summarize(mn=mean(FAGE,na.rm=T)) %>% filter(!is.na(mn))
  
  d$FAGE[d$STUDYID=="ki1119695-PROBIT"] <- d$FAGE[d$STUDYID=="ki1119695-PROBIT"]/365.25
  
  d %>% group_by(STUDYID) %>% summarize(mn=mean(FAGE,na.rm=T)) %>% filter(!is.na(mn))
  
  
  #------------------------------------------------
  # 17) Education
  #------------------------------------------------
  d$FEDUCYRS <- as.numeric(d$FEDUCYRS)
  table(d$FEDUCYRS)
  
  ################################################# 
  # Household characteristics
  #################################################
  
  #------------------------------------------------
  # 18) Asset-based wealth index
  #------------------------------------------------
  
  #Merge in HHwealth -asset based PCA:
  load("U:/data/allGHAPstudies-HHwealth.Rdata")
  colnames(pca)
  d<-left_join(d, pca, by=c("STUDYID", "COUNTRY" ,"SUBJID"))
  table(d$HHwealth_quart)
  
  
  #merge in SES variable
  table(d$SES, d$STUDYID)
  table(d$HHwealth_quart, d$STUDYID)
  d$SES[d$SES=="Upper-middle"]<-"Upper"
  d$HHwealth_quart[is.na(d$HHwealth_quart) & d$SES=="Low"] <- "Wealth Q1"
  d$HHwealth_quart[is.na(d$HHwealth_quart) & d$SES=="Lower-middle"] <- "Wealth Q2"
  d$HHwealth_quart[is.na(d$HHwealth_quart) & d$SES=="Middle"] <- "Wealth Q3"
  d$HHwealth_quart[is.na(d$HHwealth_quart) & d$SES=="Upper"] <- "Wealth Q4"
  table(d$HHwealth_quart)
  
  #income
  d$INCTOT <- as.numeric(d$INCTOT)
  summary(d$INCTOT)
  
  
  #------------------------------------------------
  # 19) Number of rooms
  #------------------------------------------------
  
  d$NROOMS <- as.numeric(d$NROOMS)
  table(d$NROOMS)
  quantile(d$NROOMS, na.rm=T)
  d$nroom<-NA
  d$nroom[d$NROOMS<2] <- "1"
  d$nroom[d$NROOMS==2] <- "2"
  d$nroom[d$NROOMS==3] <- "3"
  d$nroom[d$NROOMS>3] <- "4+"
  d$nroom <- as.factor(d$nroom)
  table(d$nroom)
  
  #------------------------------------------------
  # 20) Total number of persons
  #------------------------------------------------
  
  d$NCOMP <- as.numeric(d$NCOMP)
  table(d$NCOMP)
  
  quantile(d$NCOMP, na.rm=T)
  d$ncomp<-NA
  d$ncomp[d$NCOMP<6] <- "5 or less"
  d$ncomp[d$NCOMP>5 & d$NCOMP<9] <- "6-8"
  d$ncomp[d$NCOMP>8 & d$NCOMP<13] <- "9-12"
  d$ncomp[d$NCOMP>12] <- "13+"
  d$ncomp <- as.factor(d$ncomp)
  table(d$ncomp)
  
  #------------------------------------------------
  # 21) Number of children <5y
  #------------------------------------------------
  
  d$NCHLDLT5 <- as.numeric(d$NCHLDLT5)
  table(d$NCHLDLT5)
  #change to categorical
  d$nchild5<-NA
  d$nchild5[d$NCHLDLT5==0] <- "0"
  d$nchild5[d$NCHLDLT5==1] <- "1"
  d$nchild5[d$NCHLDLT5>2] <- "2+"
  d$nchild5 <- as.factor(d$nchild5)
  table(d$nchild5)
  
  #------------------------------------------------
  # 22) Animals (owns chickens, owns cow/buffalo)
  #------------------------------------------------
  
  #chicken
  d$CHICKEN[d$CHICKEN=="Y"]<-"1"
  d$CHICKEN<-as.numeric(as.character(d$CHICKEN))
  table(d$CHICKEN)
  d$chicken<-ifelse(d$CHICKEN==0,0,1)
  d$chicken[d$CHICKEN==""]<-NA
  table(d$chicken)
  
  #cattle
  d$CATTLE[is.na(d$CATTLE)]<-d$COW[is.na(d$CATTLE)]
  table(d$CATTLE)
  
  d$CATTLE<-as.character(d$CATTLE)
  d$CATTLE[d$CATTLE=="Yes"]<-"1"
  d$CATTLE[d$CATTLE=="No"]<-"0"
  d$CATTLE[d$CATTLE==" 1"]<-"1"
  d$CATTLE[d$CATTLE==" 0"]<-"0"
  d$CATTLE[d$CATTLE==""]<-"NA"
  d$CATTLE<-as.numeric(as.character(d$CATTLE))
  
  d$cow<-ifelse(d$CATTLE>0,1,0)
  d$cow[is.na(d$CATTLE)]<-NA
  table(d$cow)
  
  
  #------------------------------------------------
  # 23) Improved floor (vs. unimproved)
  #------------------------------------------------
  
  table(d$FLOOR)
  d$improved.floor<-ifelse(d$FLOOR=="Clay" | d$FLOOR=="Clay and Dung and Clay" | d$FLOOR=="Clay and Sand" | d$FLOOR=="Dirt as floor" | d$FLOOR=="Dung" |
                             d$FLOOR=="Dung and Clay" | d$FLOOR=="Dung and Clay and Sand" | d$FLOOR=="earth" | d$FLOOR=="Earth" | d$FLOOR=="Earth or bamboo" |
                             d$FLOOR=="Earth/Sand" | d$FLOOR=="ground" | d$FLOOR=="Ground" | d$FLOOR=="HALF SEND AND HALF CEMMENTED." | d$FLOOR=="HOUSE BEING BUILD" |
                             d$FLOOR=="Katcha" | d$FLOOR=="Katcha/Mud" | d$FLOOR=="MUD,GRAVEL" | d$FLOOR=="ONE ROOM HAS CEMENT AND ONE SAND (FLOOR)" |
                             d$FLOOR=="Palm/Bamboo" | d$FLOOR=="Thatch, grass, sticks, branches" | d$FLOOR=="cane" | d$FLOOR=="0"
                           ,0,1)
  
  d$improved.floor[d$FLOOR == "F"  | d$FLOOR == "Does not know" | d$FLOOR == "No walls/ Fence" | d$FLOOR=="Other" | d$FLOOR=="other" | d$FLOOR==""]<-NA
  
  
  #------------------------------------------------
  # 24) Sanitation facility (JMP definitions)
  #------------------------------------------------
  #Replace sanitation with sanitation arm from WASH Benefits studies
  table(d$ARM)
  
  d$SANITATN[d$shortid=="wsb" | d$shortid=="wsk"]<-"No sanitation intervention"
  d$SANITATN[d$ARM=="WSH" | d$ARM=="Sanitation" | d$ARM=="Nutrition + WSH"]<-"Sanitation intervention"
  table(d$SANITATN)
  
  
  d$improved.sanitation<-ifelse(
    d$SANITATN=="Drain connected inside house" |
      d$SANITATN=="Flush to piped sewer system" |
      d$SANITATN=="Flush to septic tank." |
      d$SANITATN=="Flush toilet" |
      d$SANITATN=="Inside latrine with drainage" |
      d$SANITATN=="Latrine with flush system" |
      d$SANITATN=="Own flush latrine" |
      d$SANITATN=="Septic tank or toilet" |
      d$SANITATN=="Latrine Septic tank/ Modern toilet" |
      d$SANITATN=="Water Closet" |
      d$SANITATN=="1" |    d$SANITATN==" 1" |
      d$SANITATN=="concrete slb" |
      d$SANITATN=="concrete slb|potty" |
      d$SANITATN=="concrete slb|waterseal" |
      d$SANITATN=="concrete slb|waterseal|potty" |
      d$SANITATN=="Own latrine|concrete slb" |
      d$SANITATN=="Own latrine|concrete slb|potty" |
      d$SANITATN=="Own latrine|concrete slb|waterseal" |
      d$SANITATN=="Own latrine|concrete slb|waterseal|potty" |
      d$SANITATN=="Own latrine|waterseal" |
      d$SANITATN=="Pacca latrine (water seal)" |
      d$SANITATN=="Sanitary" |
      d$SANITATN=="VIP Latrine w/ Water Seal" |
      d$SANITATN=="Water-sealed or slab latrine" |
      d$SANITATN=="Water sealed/slab" |
      d$SANITATN=="waterseal" |
      d$SANITATN=="Y"|
      d$SANITATN=="Sanitation intervention"|
      d$SANITATN=="toilet" |
      d$SANITATN=="TOILET" |
      d$SANITATN=="Vent. Impr. pit latrine" |
      d$SANITATN=="Vent.impr.pit latrine" |
      d$SANITATN=="VIP Latrine" |
      d$SANITATN=="Pour Flush Toilet",
    1,0)
  
  d$improved.sanitation[d$SANITATN=="" | d$SANITATN=="9" | d$SANITATN==" 9" | d$SANITATN=="Other" | d$SANITATN=="other"]<-NA
  
  #Merge in improved latrine indicator from WASH Benefits
  #"IMPRLAT"
  
  
  #------------------------------------------------
  # 25) Safe water source (JMP definitions)
  #------------------------------------------------
  
  #Replace water source with water arm from WASH Benefits studies
  table(d$ARM)
  d$H2OSRC[d$shortid=="wsb" | d$shortid=="wsk"]<-"No water intervention"
  d$H2OSRC[d$ARM=="WSH" | d$ARM=="Water" | d$ARM=="Nutrition + WSH"]<-"Safe water intervention"
  table(d$H2OSRC)
  
  
  table(d$H2OSRC)
  d$safe.water<-ifelse(
    d$H2OSRC %in%
      c("Safe water intervention",
        "Deep Bore Hole and Public Tap",
        "Public Handpump and Public Tap",
        "Private Handpump and Public Tap", 
        "Pipe and Public Tap",
        "Deep Bore Hole and Public Tap",
        "Private Well",
        "Public Handpump and Public Tap",
        "Public Tap",
        "Cover well in neigbours",
        "Bore hole tap water into neighbor's comp",
        "Covered well from neighbor's comp", 
        "Tap water into neighbor's comp",
        "Tap water from neighbor's comp",   
        "Hand pump",                 
        "Public tap",
        "Taken from mosque next door through pipe",
        "Cover well in neigbours", 
        "Bore hole tap water into neighbor's comp",
        "Covered well from neighbor's comp",
        "Tap water into neighbor's comp",
        "Tap water from neighbor's comp",
        "Covered well into neighbours camp",
        "Tap in neighbor's compound",
        "Pipe into neighbor's yard",
        "Pipe into neighbor hosue", 
        "Tap in neighbours compound",
        "Piped into relative house",
        "Borrowed from neighbour through pipe",
        "From waterpump",                         
        "Piped from neighbors",   
        "Land lord's tap",  
        "Rainwater",
        "Piped water at the neighbors",           
        "Bore hole at neighbors",  
        "Protected spring",
        "Tap into a neighbors underground", 
        "Tap service office",
        "Covered private well",
        "Tap in house",                            "Connection in Home for 24hrs",           
        "Other Source",                           
        "Private tap",                             "Community tap",                          
        "Private tubewell/ hand pump",             "Community tubewell",                     
        "Piped into dwelling from treated faclty",
        "Piped into dwelling potable",             "Piped to compnd, yard or plot-trt faclty",
        "Piped to compnd, yard or plot-potable",                         "Piped to neighbor -potable",             
        "public tap-treatment faclty",             "dug well- protected",                    
        "municipal network",                                                       
        "Deep tube well",                         
        "Piped into yard",                        
        "Bore hole",                               "Covered public well",                    
        "Piped into neighbors compound",          
        "Piped into house",                                 
        "Piped into nearby compound",              "Covered well in house or yard")
    ,1,0)
  table(d$safe.water)
  
  
  d$safe.water[d$H2OSRC=="Other" | d$H2OSRC=="other" | d$H2OSRC=="Other Source" | d$H2OSRC==""
               | d$H2OSRC=="Self-made waterline is near: N/A" | d$H2OSRC=="Self-made waterline is near: Drain" | 
                 d$H2OSRC=="Self-made waterline is near: Latrine" | d$H2OSRC=="Drinking water" ]<-NA
  table(d$safe.water)
  
  #------------------------------------------------
  # 26) Treats drinking water
  #------------------------------------------------
  
  #Replace water treatment with water arm from WASH Benefits studies
  table(d$ARM)
  d$H2OTRTP[d$STUDYID=="ki1000110-WASH-Bangladesh" | d$STUDYID=="ki1000111-WASH-Kenya"]<-"No water intervention"
  d$H2OTRTP[d$ARM=="WSH" | d$ARM=="Water" | d$ARM=="Nutrition + WSH"]<-"Safe water intervention"
  table(d$H2OTRTP)
  
  d$treat.water<-ifelse(
    d$H2OTRTP=="Strain Through Cloth" | 
      d$H2OTRTP=="No water intervention" |
      d$H2OTRTP=="Soda" | 
      d$H2OTRTP=="None" | 
      d$H2OTRTP=="Leave water in sun" | 
      d$H2OTRTP=="Let It Stand" | 
      d$H2OTRTP=="Filter through a cloth" | 
      d$H2OTRTP=="Do not treat" | 
      d$H2OTRTP=="Not Boiled" | 
      d$H2OTRTP=="No treatment" | 
      d$H2OTRTP=="Alum" | 
      d$H2OTRTP=="Own Arrangement by plastic pipe" | 
      d$H2OTRTP=="Municipality supply" | 
      d$H2OTRTP=="Alum" 
    ,0,1)
  
  table(d$treat.water)
  
  d$treat.water[d$H2OTRTP=="does not know" | d$H2OTRTP=="unknown or not applicable" | d$H2OTRTP=="Other"  | d$H2OTRTP=="Others"| d$H2OTRTP==""]<-NA
  
  
  
  #Merge in the unsafe water variable
  table(d$H2OUNTRT)
  
  
  #------------------------------------------------
  # 27) Soap present in the home
  #------------------------------------------------
  
  table(d$SOAP)
  
  #Add in uptake indicator from WASH B
  #"WATSOAP"
  
  #------------------------------------------------
  # 28) Cooking fuel usage
  #------------------------------------------------
  
  table(d$COOKFUEL)
  
  d$cleancook <- ifelse(d$COOKFUEL=="Electricity"|
                          d$COOKFUEL=="Gas"|
                          d$COOKFUEL=="Gas cylinder stove"|
                          d$COOKFUEL=="gas stove"|
                          d$COOKFUEL=="Gas Stove"|
                          d$COOKFUEL=="Kerosene oil stove"|
                          d$COOKFUEL=="kerosene stove"|
                          d$COOKFUEL=="more than one",1,0)
  d$cleancook[is.na(d$COOKFUEL)] <- NA
  d$cleancook[d$COOKFUEL=="Other" | d$COOKFUEL=="Garments products"] <- NA
  table(d$cleancook)
  
  
  #------------------------------------------------
  # 28.5) Food security
  #------------------------------------------------
  
  table(d$FOODDFCT)
  
  
  
  #################################################
  # Time-varying characteristics
  #################################################
  
  #------------------------------------------------
  # 29) Current diarrhea
  #------------------------------------------------
  d$DIARFL[d$DIARFL==" 0"]<-"0"
  d$DIARFL[d$DIARFL==" 1"]<-"1"
  d$DIARFL<-as.numeric(d$DIARFL)
  
  #------------------------------------------------
  # 30) Cumulative days of diarrhea
  #------------------------------------------------
  #create harmonized diarrhea variable -add in maled
  table(d$STUDYID, is.na(d$PCTDIAR))
  d$PCTDIAR[d$STUDYID=="MAL-ED"]<-as.numeric(as.character(d$DIARDAYS[d$STUDYID=="MAL-ED"]))/30*100
  
  #Other variables
  #"SUMEP","SUMDIAR",
  #"SUMDAYS",
  
  #"DIARFL","SUMDIAR",  
  #"DIARDAYS", "CSUMDIAR",
  #"DIAREPS" , "DIARBFEN" ,
  #"DIARRHOEA","DIARRHOEA_NEONATAL"
  
  #------------------------------------------------
  # 31) Breastfeeding duration
  #------------------------------------------------
  
  #"DURBRST"
  table(d$DURBRST)
  d$DURBRST <- as.numeric(d$DURBRST)
  
  #------------------------------------------------
  # 32) Season (Month of measurement)
  #------------------------------------------------
  
  #"BRTHWEEK"
  
  d$month <- floor(d$BRTHWEEK/4 + d$AGEDAYS/30.25)
  table(d$month)
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  d$month[d$month>12 & !is.na(d$month)] <- d$month[d$month>12 & !is.na(d$month)]-12 
  table(d$month)
  
  
  #################################################
  # Intervention ARM
  #################################################
  d$ARM<-as.character(d$ARM)
  
  #Replace nondescript intervention names
  d$ARM[d$STUDYID=="ki1000125-AgaKhanUniv" & d$ARM=="Intervention"] <-"Maternal Education"
  d$ARM[d$STUDYID=="ki1000304b-SAS-CompFeed" & d$ARM=="Intervention"] <-"Comp. feeding education"
  
  d$ARM[d$STUDYID=="ki1148112-iLiNS-DYAD-M" & d$ARM=="Multiple micronutrient supplementation"] <-"Maternal Multiple micronutrient supplementation"
  d$ARM[d$STUDYID=="ki1148112-iLiNS-DYAD-M" & d$ARM=="Lipid-based nutrient supplementation"] <-"Maternal and child LNS"
  d$ARM[d$STUDYID=="kiGH5241-JiVitA-3" & d$ARM=="Multiple Micronutrients"] <-"Maternal Multiple Micronutrients"
  
  
  d$tr <- NA
  #Zinc: Supplement zinc, or food fortified with zinc. 
  d$tr[d$ARM=="Therapeutic Zinc: 20 mg/day for 10 days" | d$ARM=="3 mg zinc, no copper" | d$ARM=="10 mg zinc, no copper" |  d$ARM=="10 mg zinc, with copper" |  d$ARM=="3 mg zinc, no copper" | 
         d$ARM=="Intermittent Zinc: 10 mg/d for 10 days" | d$ARM=="Preventive Zinc: 7 mg/day" | d$ARM=="Zinc Alone" | d$ARM=="7 mg zinc, no copper" | d$ARM=="Zinc" |
         d$ARM=="b.LNS-Zn5" | d$ARM=="c.LNS-Zn10" | d$ARM=="d.LNS-TabZn5"  ] <- "Zinc"
  #Need to write code to set a.LNS-Zn0 as control arm in LNS zinc trial
  
  #Lipid based nutrient supplements (LNS): Nutritional supplements that deliver to the child energy, protein, and fatty acids, mostly from lipids, as well as micronutrients.
  d$tr[d$ARM=="a.LNS-Zn0" | d$ARM=="LNS-20gNoM" | d$ARM=="LNS-20gM" | d$ARM=="LNS-10gM" | d$ARM=="LNS-40gM" | d$ARM=="LNS-40gNoM" | 
         d$ARM=="Nutrition" | d$ARM=="Nutrition + WSH"   |  d$ARM=="Plumpy Doz" | d$ARM=="Lipid-based nutrient supplementation" | d$ARM=="Maternal and child LNS" |  
         d$ARM=="Milk FS" | d$ARM=="Soy FS" | d$ARM=="Chickpea" | d$ARM=="Rice Lentil"] <- "LNS"
  
  #Maternal interventions: Prenatal interventions delivered to the mother such as maternal multiple micronutrient supplementation, LNS, or education interventions
  d$tr[d$ARM=="Maternal Education" | d$ARM=="Maternal and child LNS" | 
         d$ARM=="Maternal Multiple micronutrient supplementation" | d$ARM=="Maternal Multiple Micronutrients" ] <- "Mat"
  
  #Other interventions: Interventions that do not belong to any of the above categories. No pooled effects of other interventions will be calculated; they will be analyzed following the methods laid out in this plan but interpreted based on study-specific intervention design. Other interventions include multiple and single (other than zinc) micronutrient supplementation, complementary feeding, and child nutritional counseling.
  d$tr[d$ARM=="WSH" | d$ARM=="Water" | d$ARM=="Handwashing" | d$ARM=="Sanitation" | d$ARM=="Education" | d$ARM=="Visitation" | d$ARM=="Nutritional counselling" | d$ARM=="Vitamin D"  | 
         d$ARM=="50,000 IU nippled + 400,000 IU Oval" | d$ARM=="Placebo nippled + 400,000 IU Oval" | d$ARM=="50,000 IU nippled + Placebo Oval" | d$ARM=="BSC" | d$ARM=="Comp. feeding education" | d$ARM=="Nutritional counselling" | d$ARM=="Visitation"] <- "Other"
  
  d$tr[ d$ARM=="Multivitamin Alone" | d$ARM=="Zinc + Multivitamin" |  d$ARM=="MNT + WPC" | d$ARM=="MNT + BSC"] <- "Other" #MMN interventions
  
  d$tr[d$ARM=="Food supplementation"  | d$ARM=="WSB++" ] <- "Other" #Complementary feeding interventions
  
  #Control arms: No intervention, the standard of care in the study region (i.e. maternal iron and folate supplementation), or the intervention arm presented as the control in the trial design. If the intervention is factorial in design (i.e. iLiNS Zinc), the control arm will be chosen to isolate the effect of the intervention type (so in iLiNS Zinc, arms receiving LNS and Zinc will be compared to the LNS-only arm rather than the control arm).
  d$tr[d$ARM=="Control" | d$ARM=="Control (no Zinc)" | d$ARM=="Standard(Control)" | d$ARM=="No intervention" | d$ARM=="Placebo" | d$ARM=="Passive Control" | d$ARM=="no zinc, no copper" 
       | d$ARM=="Iron and folic acid supplementation" | d$ARM=="e.Control" | d$ARM=="Likuni Phala" | d$ARM=="WPC" | d$ARM== "CFC"  | d$ARM=="Placebo nippled + Placebo Oval" | d$ARM=="Iron Folic Acid"] <- "C"
  
  d$tr <- factor(d$tr)
  
  
  
  
  #################################################
  # Others not yet registered
  #################################################
  
  
  #Bednet
  #d$bednet<-ifelse(d$BEDNET==" 0",0,1)
  #d$bednet[is.na(d$BEDNET)]<-NA
  
  
  
  
  
  
  
  #Create a region of the world variable
  unique(d$COUNTRY)
  d$region<-"Missing"
  d$region[d$COUNTRY=="BURKINA FASO" | d$COUNTRY=="GAMBIA" | d$COUNTRY=="GUINEA-BISSAU" 
           | d$COUNTRY=="KENYA" | d$COUNTRY=="MALAWI" | d$COUNTRY=="MALI" | d$COUNTRY=="MOZAMBIQUE"
           | d$COUNTRY=="SOUTH AFRICA" | d$COUNTRY=="TANZANIA, UNITED REPUBLIC OF" | d$COUNTRY=="ZIMBABWE"]<-"SSA"
  d$region[d$COUNTRY=="BANGLADESH" | d$COUNTRY=="INDIA" | d$COUNTRY=="NEPAL" | d$COUNTRY=="PAKISTAN" ]<-"south Asia"
  d$region[d$COUNTRY=="CHINA" | d$COUNTRY=="PHILIPPINES" | d$COUNTRY=="NEPAL" | d$COUNTRY=="PAKISTAN" ]<-"east Asia"
  d$region[d$COUNTRY=="BRAZIL" | d$COUNTRY=="ECUADOR" | d$COUNTRY=="GUATEMALA" | d$COUNTRY=="PERU" ]<-"Latin America"
  d$region[d$COUNTRY=="BELARUS"]<-"Europe"
  table(d$region)
  d$region<-as.factor(d$region)
  d$COUNTRY<-as.factor(d$COUNTRY)
  
  
  
  
  
  
  
  d<- d[,which(colnames(d) %in% c(
    "STUDYID",
    "COUNTRY",
    "SUBJID",
    "AGEDAYS",
    "HAZ",
    "stunt",
    "sevstunt",
    "WHZ",   "wast","sevwast","period_length","sevwast_inc",
    "wast_inc","sevwast_rec","wast_rec","wast_risk",
    "sevwast_risk","wast_rec_risk","sevwast_rec_risk","wasting_episode",
    "born_wast_inc","episode_ID","incident_age","maxage",
    "duration","wasting_duration","born_sevwast_inc","sevwasting_episode", 
    "sev_episode_ID","sevwasting_duration", "wast_rec90d","wast_rec60d",
    "wast_rec30d","sevwast_inc90d","sevwast_inc60d","sevwast_inc30d",  
    "SEX",
    "GAGEBRTH",
    "BIRTHLEN",
    "BIRTHWT",
    "birthorder",
    "birthmonth",
    "homedelivery",
    "vagbirth",
    "MHTCM",
    "MWTKG", 
    "MBMI", 
    "MAGE", 
    "MEDUCYRS",
    "single",
    "breastfeeding",
    "FHTCM", 
    "FAGE",
    "FEDUCYRS",
    "nroom",
    "ncomp",
    "nchild5",
    "chicken",
    "cow",
    "improved.floor",
    "improved.sanitation",
    "safe.water",
    "treat.water",
    "SOAP",
    "cleancook",
    "DIARFL",
    "DURBRST", 
    "month",
    "tr", "ARM",
    "enrolstunt",
    "HHwealth_quart",
    "region"))]
  
  
  # Impute missing data with cohort-specific medians and make flagged indicators for missing values
  #library(ck37r)
  
  d$birthorder <- factor(d$birthorder)
  d$SOAP <- factor(d$SOAP)
  d$STUDYID <- factor(d$STUDYID)
  
  d<-as.data.frame(d)
  
  #Convert outcomes to numeric
  outcomes<-c("stunt",
              "sevstunt",
              "WHZ",   "wast","sevwast","period_length","sevwast_inc",
              "wast_inc","sevwast_rec","wast_rec","wast_risk",
              "sevwast_risk","wast_rec_risk","sevwast_rec_risk","wasting_episode",
              "born_wast_inc","episode_ID","incident_age","maxage",
              "duration","wasting_duration","born_sevwast_inc","sevwasting_episode",
              "sev_episode_ID","sevwasting_duration", "wast_rec90d","wast_rec60d",
              "wast_rec30d","sevwast_inc90d","sevwast_inc60d","sevwast_inc30d")
  
  for(i in outcomes){
    if(i %in% colnames(d)){
      d[,which(colnames(d) %in% i)]<-as.numeric(as.character(d[,which(colnames(d) %in% i)]))
    }
  }
  
  
  return(d)
  
}






