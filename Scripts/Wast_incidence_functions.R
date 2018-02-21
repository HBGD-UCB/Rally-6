

rm(list=ls())
library(tidyverse)
library(zoo)



#Mock data
# set.seed(12345)
# AGEDAYS<-c(1,31,70, 100, 200, 210,220,270,300,
#            5,20,50,100,130,
#            7,14,21,28,45)
# AGEDAYS<-rep(sample.int(720, 60),50)
# SUBJID<-c(1,1,1,1,1, 1,1,1,1,
#           2,2,2,2,2,
#           3,3,3,3,3)
# SUBJID<-rep(1:50, each = 60)
# WHZ<-c(0,-3.5,0,-3.5,0,0,0,-3.5,0,
#        -2.5,-2.5,-2.5,-2.5,-2.5,
#        -2.5,0,-2.5,0,0)
# WHZ<-rnorm(3000,-2,1)
# d<-data.frame(SUBJID,AGEDAYS,WHZ)
# d<-d %>% arrange(SUBJID,AGEDAYS)
# 
# agecats=c(6*30, 12*30, 18*30, 24*30)
#   d$agecat <- as.factor(findInterval(d$AGEDAYS, agecats, rightmost.closed=F))
# table(d$agecat)


# strat=T
# agecats=c(6*30, 12*30, 18*30, 24*30)
# agecat_rownames=NULL
# washout=60
# 
# test<- WastIncCalc(d,washout=60)
# test2<-WastIncTable(test,
#                     strat=T,
# agecats=c(6*30, 12*30, 18*30, 24*30),
# agecat_rownames=NULL)



#----------------------------------------------
#create functions for rolling sum windows 
#----------------------------------------------

roll_sum_fun <-  function(v, len){ sapply(1:(length(v)),function(x){sum(v[(x+1):(x+len+1)], na.rm=T)})}
lag_sum_fun <-  function(v, len){ sapply(1:(length(v)),function(x){ifelse((x-1-len)<0,
                                                                          sum(v[0:(x-1)], na.rm=T),
                                                                          sum(v[(x-1-len):(x-1)], na.rm=T))})}
#function to always round 0.5 up
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


#----------------------------------------------
# Mean and 95% CI function 
#----------------------------------------------
mean95CI <- function(Y, id=rep(1:length(Y)), persontime=NULL, proportion=F, percent=F, count=F){
  
  if(proportion==F){
    if(count==T){
      IR.CI <- pois.exact(Y, pt = persontime, conf.level = 0.95)[3:5] 
      mean_ci <- data.frame(N=Y, Mean=IR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=IR.CI[2] ,  Upper.95.CI=IR.CI[3] )
      colnames(mean_ci) <- c("N","Mean","SD","Robust SE", "Lower 95%CI", "Upper 95%CI") 
    }else{
      if(!is.na(mean(Y[complete.cases(Y)]))){
        mudat <- data.frame(id = id, Y = Y)
        mudat <- mudat[complete.cases(mudat), ]
        n.sub <- dim(mudat)[1]
        fit <- glm(Y ~ 1, family = gaussian, data = mudat)
        vcovCL <- sandwichSE(mudat, fm = fit, cluster = mudat$id)
        rfit <- coeftest(fit, vcovCL)
        lb <- rfit[1, 1] - 1.96 * rfit[1, 2]
        ub <- rfit[1, 1] + 1.96 * rfit[1, 2]
        mean_ci <- matrix(c(n.sub, rfit[1, 1], sd(mudat$Y), rfit[1, 
                                                                 2], lb, ub), nrow = 1, ncol = 6)
        colnames(mean_ci) <- c("N", "Mean", "SD", "Robust SE", "Lower 95%CI", 
                               "Upper 95%CI")
        
      }else{
        mean_ci <- data.frame(N=NA, Mean=NA, SD=NA, `Robust SE`=NA, `Lower 95%CI`=NA, `Upper 95%CI`=NA)
        colnames(mean_ci) <- c("N", "Mean", "SD", "Robust SE", "Lower 95%CI", "Upper 95%CI")  
      }
    }
  }else{
    
    require(binom)
    # Find the number of obs
    n = length(Y[!is.na(Y)])
    if(percent==T){
      CR.res<-binom.confint(sum(Y/100, na.rm = T), n, method="exact")
    }else{
      CR.res<-binom.confint(sum(Y, na.rm = T), n, method="exact")
    }
    mean_ci <- data.frame(N=n, Mean=CR.res[4], SD=NA, `Robust SE`=NA, `Lower 95%CI`=CR.res[5], `Upper 95%CI`=CR.res[6])
    colnames(mean_ci) <- c("N", "Mean", "SD", "Robust SE", "Lower 95%CI", "Upper 95%CI")  
  }
  return(mean_ci)
}



#Function to calculate the robust SEs
sandwichSE <- function (dat, fm, cluster) 
{
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  if (is.factor(cluster)) {
    cluster <- droplevels(cluster)
  }
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj <- apply(estfun(fm), 2, function(x) tapply(x, cluster, 
                                                sum))
  vcovCL <- dfc * sandwich(fm, meat = crossprod(uj)/N)
  return(vcovCL)
}


#----------------------------------------------
#create function to calc unstrat and age stat incidence
#----------------------------------------------
WastIncCalc<-function(d, washout=60, dropBornWasted=F){
  require(tidyverse)
  require(zoo)  
  
  #Filter out extreme or missing whz values
  d <- d %>%  filter(!is.na(WHZ)) %>%
    filter(WHZ > (-5) & WHZ < 5)
  
  #Remove duplicate ages
  ndropped <- nrow(d[duplicated(cbind(d$SUBJID, d$AGEDAYS)), ])
  d <- d[!duplicated(cbind(d$SUBJID, d$AGEDAYS)), ]
  if(ndropped>0) cat("\n-----------------------------------\n",ndropped," observations dropped due to duplicate ages\n_----------------------------------\n")
  
  
  #Create visit variable
  #(This replaces the visit variables, but some studies have missing 
  #visit data)
  d <- d %>% group_by(SUBJID) %>% mutate(VISITNUM = rank(AGEDAYS))
  
  
  #Extract required columns and save others to merge back in later
  othercolumns <- d %>% subset(., select= -c(WHZ, VISITNUM)) 
  d <- d %>% subset(., select= c(SUBJID, WHZ, AGEDAYS, VISITNUM)) 
  
  
  #generate wasting and severe wasting indicators
  d$wast= ifelse(d$WHZ < (-2),1,0)
  d$sevwast= ifelse(d$WHZ < (-3),1,0)
  
  # #Generate variables for length of period in days between prior observation and current observations
  # #and the next observations and current observations. Also generate variables for if child changed from
  # #not wasted to wasted (or severe wasted) between the last observation and the current observation.
  d <- d %>%
    arrange(SUBJID, AGEDAYS) %>%
    group_by(SUBJID) %>%
    mutate(
      agelag=lag(AGEDAYS),
      wastlag=lag(wast),
      sevwastlag=lag(sevwast),
      midpoint_age = AGEDAYS - (AGEDAYS - agelag)/2,
      wastchange = wast - lag(wast),
      sevwastchange = sevwast - lag(sevwast)
    ) %>%
    as.data.frame()
  d$agelag[is.na(d$agelag)] <- 0
  d$wastlag[is.na(d$wastlag)] <- 0
  d$sevwastlag[is.na(d$sevwastlag)] <- 0
  d$wastchange[is.na(d$wastchange)] <- d$wast[is.na(d$wastchange)]
  d$sevwastchange[is.na(d$sevwastchange)] <- d$sevwast[is.na(d$sevwastchange)]
  d$midpoint_age[is.na(d$midpoint_age)] <- d$AGEDAYS[is.na(d$midpoint_age)]/2
  
  #Length of each observation period
  d <- d %>% group_by(SUBJID) %>% 
    mutate(
      next_midpoint = lead(midpoint_age)
    )
  
  #Assume 30 day period after final measurement (so midpoint will be 15 days after final measure)
  d$next_midpoint[is.na(d$next_midpoint)] <- d$AGEDAYS[is.na(d$next_midpoint)] + 15
  d$period_length <- (d$next_midpoint - d$midpoint_age)
  
  
  
  N <- nrow(d)
  d$washout_period_lead <- d$washout_period_lag <- rep(T, N)
  d$future_sevwast <- d$future_wast <-  d$past_sevwast <- d$past_wast <- rep(0, N)
  
  for(i in 1:washout){
    
    d <- d %>% group_by(SUBJID) %>%
      mutate(
        wast_lag_i  = lag(wast, i),
        sevwast_lag_i  = lag(sevwast, i),
        days_lag_i = abs(lag(midpoint_age, i) - midpoint_age),
        wast_lead_i  = lead(wast, i),
        sevwast_lead_i  = lead(sevwast, i),
        days_lead_i = abs(lead(midpoint_age, i) - midpoint_age)
      )
    
    d$washout_period_lag[d$days_lag_i > washout] <- F
    d$washout_period_lead[d$days_lead_i > washout] <- F
    
    d$past_wast[d$wast_lag_i==1 & d$washout_period_lag==T] <- 1
    d$past_sevwast[d$sevwast_lag_i==1 & d$washout_period_lag==T] <- 1              
    d$future_wast[d$wast_lead_i==1 & d$washout_period_lead==T] <- 1
    d$future_sevwast[d$sevwast_lead_i==1 & d$washout_period_lead==T] <- 1
    
    #Stop for loop if all current leading and lagging observations are beyond washout period
    if(min(d$days_lead_i, na.rm=T) & min(d$days_lag_i, na.rm=T) > washout) break
  }
  
  d <- d %>% 
    subset(., select= -c(washout_period_lag, washout_period_lead, 
                         wast_lag_i, sevwast_lag_i, days_lag_i, 
                         wast_lead_i, sevwast_lead_i, days_lead_i)) %>% 
    ungroup() %>% as.data.frame()
  
  head(d,30)
  
  #---------------------------------------------------------
  #Calculate wasting and wasting recovery incidence and risk
  #---------------------------------------------------------
  d$wast_rec_inc <- d$sevwast_rec_inc <- d$wast_inc <- d$sevwast_inc <- rep(0, N)
  d$sevwast_rec_risk <- d$wast_rec_risk <- d$sevwast_risk <- d$wast_risk <-  rep(0, N)
  
  
  d$wast_inc[d$wastchange==1 & d$past_wast==0] <- 1 #Wasting incidence if at risk of wasting and change in status between prior and current observation
  d$wast_rec_inc[d$wastchange== -1 & d$future_wast==0] <- 1 #Recovery from wasting if status change to not wasted and no new wasting in the future washout period
  
  d$wast_inc[d$wastchange!=1 | d$past_wast==1] <- 0 
  d$wast_rec_inc[d$wastchange!= -1 | d$future_wast==1] <- 0  
  
  #Remove incidences of wasting if there has not been recovery from prior wasting episode
  #Is there a cleaner way of preventing recording the incidences earlier?
  table(d$wast_inc)
  
  d <- d %>% group_by(SUBJID) %>% 
    mutate(sum_wast_inc=cumsum(wast_inc),
           sum_wast_rec=cumsum(wast_rec_inc))
  for(i in 1:nrow(d)){
    if(d$wast_inc[i]==1 & (d$sum_wast_inc[i]-d$sum_wast_rec[i] > 1)){
      d$wast_inc[i] <- 0 
      d <- d %>% group_by(SUBJID) %>% 
        mutate(sum_wast_inc=cumsum(wast_inc),
               sum_wast_rec=cumsum(wast_rec_inc))        
    }
    if(d$wast_rec_inc[i]==1 & (d$sum_wast_inc[i]-d$sum_wast_rec[i] < 0)){
      d$wast_rec_inc[i] <- 0 
      d <- d %>% group_by(SUBJID) %>% 
        mutate(sum_wast_inc=cumsum(wast_inc),
               sum_wast_rec=cumsum(wast_rec_inc))  
    }
  }
  table(d$wast_inc)
  
  
  d <- subset(d, select = -c(sum_wast_inc,sum_wast_rec))
  #Make sure there isn't double recovery
  # d <- d %>% group_by(SUBJID) %>% 
  #   mutate(sum_wast_inc=cumsum(wast_inc),
  #          sum_wast_rec=cumsum(wast_rec_inc))
  # table(d$sum_wast_inc- d$sum_wast_rec != 1 & d$sum_wast_inc- d$sum_wast_rec != 0)
  # d$flag <- d$sum_wast_inc- d$sum_wast_rec != 1 & d$sum_wast_inc- d$sum_wast_rec != 0
  # 
  #Indicate length of incident episodes
  d$wasting_episode <- rep(NA, N)
  d$wasting_episode[d$wast_inc==1] <- "Wasted"
  d$wasting_episode[d$wast_rec_inc==1] <- "Not Wasted"
  
  
  #Have to mark first observations as wasted or not wasted if dropBornWasted=F
  if(dropBornWasted==F){
    d <- d %>% group_by(SUBJID) %>% 
      mutate(wasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & wast==0, "Not Wasted", wasting_episode),
             wasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & wast==1, "Wasted", wasting_episode),
             born_wast_inc= 0,
             wasting_episode = na.locf(wasting_episode, fromLast=F)) %>% #Last observation carried forward 
      ungroup()
  }else{
    d <- d %>% group_by(SUBJID) %>% 
      mutate(wasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & wast==0, "Not Wasted", wasting_episode),
             wasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & wast==1, "Born Wasted", wasting_episode),
             wast_inc = ifelse(wasting_episode=="Born Wasted",0, wast_inc),
             born_wast_inc= ifelse(AGEDAYS==min(AGEDAYS) & wasting_episode=="Born Wasted",1,0),
             wasting_episode = na.locf(wasting_episode, fromLast=F)) %>% #Last observation carried forward 
      ungroup()      
  }
  
  #Indicate risk of wasting or recovery 
  d$wast_risk[(d$wasting_episode=="Not Wasted" & d$past_wast==0) | d$wast_inc==1] <- 1 
  d$wast_rec_risk[(d$wasting_episode!="Not Wasted" & d$wast_inc!=1) | d$wast_rec_inc==1] <- 1 
  
  
  
  #Calculate duration of wasting episodes
  d <- d %>%  group_by(SUBJID) %>%
    mutate(episode_ID = cumsum(born_wast_inc+wast_inc+wast_rec_inc) + 1) %>% #Create unique episode ID
    ungroup() %>% group_by(SUBJID, episode_ID) %>%
    mutate(incident_age = min(midpoint_age),
           maxage=max(AGEDAYS))
  
  
  d_episode <- d %>% 
    subset(., select=c(SUBJID, episode_ID, incident_age)) %>%
    group_by(SUBJID, episode_ID) %>%
    slice(1) %>% ungroup() %>% group_by(SUBJID) %>%
    mutate(duration=lead(incident_age)-incident_age) %>%
    subset(., select= -incident_age)
  
  d <- left_join(d, d_episode, by=c("SUBJID","episode_ID"))
  
  #Set duration of any censored episode to NA
  d <- d %>% group_by(SUBJID) %>%
    mutate(duration = ifelse(maxage==max(maxage), NA, duration)) %>% 
    ungroup()
  
  
  #Variable for duration of only wasting episodes
  d$wasting_duration <- NA
  d$wasting_duration[d$wasting_episode=="Wasted"] <- d$duration[d$wasting_episode=="Wasted"]
  
  
  
  #---------------------------------------------------------
  #Calculate severe wasting and severe wasting recovery incidence and risk
  #---------------------------------------------------------    
  
  #Mark severe wasting changes
  d$sevwast_falter <- NA
  d$sevwast_falter[d$sevwastchange==1 & d$past_sevwast==0] <- 1
  d$sevwast_falter[d$sevwastchange!=1 | d$past_sevwast==1] <- 0
  
  d$born_sevwast_inc<-0
  d$sevwasting_episode <- rep(NA, N)
  d$sevwasting_episode[d$sevwast_falter==1] <- "Severe Wasted"
  d$sevwasting_episode[d$wast_rec_inc==1] <- "Not Severe Wasted"
  
  #Have to mark first observations as wasted or not wasted if dropBornWasted=F
  if(dropBornWasted==T){
    d <- d %>% group_by(SUBJID) %>% 
      mutate(
        born_sevwast_inc= ifelse(AGEDAYS==min(AGEDAYS) & sevwast==1,1,0),
        sevwasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevwast==0, "Not Severe Wasted", sevwasting_episode),
        sevwasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevwast==1, "Born Severe Wasted", sevwasting_episode),
        sevwasting_episode = na.locf(sevwasting_episode, fromLast=F),
        sevwasting_episode_lag=lag(sevwasting_episode)) %>% #Last observation carried forward 
      ungroup()      
  }else{
    d <- d %>% group_by(SUBJID) %>% 
      mutate(
        sevwasting_episode = ifelse
        (AGEDAYS==min(AGEDAYS) & sevwast==0, "Not Severe Wasted", sevwasting_episode),
        sevwasting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevwast==1, "Severe Wasted", sevwasting_episode),
        sevwasting_episode = na.locf(sevwasting_episode, fromLast=F),
        sevwasting_episode_lag=lag(sevwasting_episode)) %>% #Last observation carried forward 
      ungroup()     
  }
  
  #Indicate incidence of severe wasting and recovery
  d$sevwasting_episode_lag[is.na(d$sevwasting_episode_lag)]<-"Not Severe Wasted"
  d$sevwast_inc <- d$sevwast_rec_inc <- 0
  d$sevwast_inc[d$sevwasting_episode=="Severe Wasted" & (d$sevwasting_episode_lag=="Not Severe Wasted")] <- 1
  d$sevwast_rec_inc[d$sevwasting_episode=="Not Severe Wasted" & (d$sevwasting_episode_lag!="Not Severe Wasted")] <- 1
  
  
  #Create unique severe wasting episode IDs
  d <- d %>%  group_by(SUBJID) %>%
    mutate(sev_episode_ID = cumsum(sevwast_inc+born_sevwast_inc+sevwast_rec_inc) + 1) %>% #Create unique episode ID
    ungroup() 
  
  #Fill in severe wasting  episodes
  # d <- d %>% group_by(SUBJID, sev_episode_ID) %>%
  #   mutate(sevwasting_episode = na.locf(sevwasting_episode, fromLast=F)) %>%
  #   ungroup()
  
  
  #Mark risk of severe wasting and severe wasting recovery 
  d$sevwast_risk[(d$sevwasting_episode!="Severe Wasted" & d$sevwasting_episode!="Born Severe Wasted" & d$past_sevwast==0) | d$sevwast_inc==1] <- 1 
  d$sevwast_rec_risk[((d$sevwasting_episode=="Severe Wasted" | d$sevwasting_episode=="Born Severe Wasted") & d$sevwast_inc!=1) | d$sevwast_rec_inc==1] <- 1 
  
  
  #Calculate duration of severe wasting episodes
  d <- d  %>% group_by(SUBJID, sev_episode_ID) %>%
    mutate(sev_incident_age = min(midpoint_age),
           sev_maxage=max(AGEDAYS))
  
  
  
  d_sev_episode <- d %>% 
    subset(., select=c(SUBJID, sev_episode_ID, sev_incident_age)) %>%
    group_by(SUBJID, sev_episode_ID) %>%
    slice(1) %>% ungroup() %>% group_by(SUBJID) %>%
    mutate(sevduration=lead(sev_incident_age)-sev_incident_age) %>%
    subset(., select= -sev_incident_age)
  
  d <- left_join(d, d_sev_episode, by=c("SUBJID","sev_episode_ID"))
  
  #Set duration of any censored episode to NA
  d <- d %>% group_by(SUBJID) %>%
    mutate(sevduration = ifelse(sev_maxage==max(sev_maxage), NA, sevduration)) %>% 
    ungroup()
  
  
  #Variable for duration of only severe wasting episodes
  d$sevwasting_duration <- NA
  if(dropBornWasted==F){
    d$sevwasting_duration[d$sevwasting_episode=="Severe Wasted"] <- d$sevduration[d$sevwasting_episode=="Severe Wasted"]
  }else{
    d$sevwasting_duration[d$sevwasting_episode=="Severe Wasted" | d$sevwasting_episode=="Born Severe Wasted"] <- d$sevduration[d$sevwasting_episode=="Severe Wasted" | d$sevwasting_episode=="Born Severe Wasted"]
  }
  
  
  #Calculate 30,60, 90 day recovery and faltering into severe wasting
  d$period_30d <- d$period_60d <- d$period_90d <- T
  d$wast_rec30d <- d$wast_rec60d <- d$wast_rec90d <- NA
  d$wast_rec30d[d$wast_inc==1]  <- d$wast_rec60d[d$wast_inc==1] <- d$wast_rec90d[d$wast_inc==1] <- 0
  d$sevwast_inc30d <- d$sevwast_inc60d <- d$sevwast_inc90d <- NA
  d$sevwast_inc30d[d$wast_inc==1] <- d$sevwast_inc60d[d$wast_inc==1] <- d$sevwast_inc90d[d$wast_inc==1] <- 0
  for(i in 1:90){
    d <- d %>% group_by(SUBJID) %>%
      mutate(
        rec_inc_lead_i = lead(wast_rec_inc, i),
        sev_inc_lead_i = lead(sevwast_inc, i),
        days_lead_i = abs(lead(midpoint_age, i) - midpoint_age)
      )
    
    d$period_30d[d$days_lead_i > 30] <- F
    d$period_60d[d$days_lead_i > 60] <- F
    d$period_90d[d$days_lead_i > 90] <- F
    
    d$wast_rec30d[d$wast_inc==1 & d$period_30d & d$rec_inc_lead_i==1] <- 1
    d$wast_rec60d[d$wast_inc==1 & d$period_60d & d$rec_inc_lead_i==1] <- 1
    d$wast_rec90d[d$wast_inc==1 & d$period_90d & d$rec_inc_lead_i==1] <- 1
    
    d$sevwast_inc30d[d$wast_inc==1 & d$period_30d & d$sev_inc_lead_i==1] <- 1
    d$sevwast_inc60d[d$wast_inc==1 & d$period_60d & d$sev_inc_lead_i==1] <- 1
    d$sevwast_inc90d[d$wast_inc==1 & d$period_90d & d$sev_inc_lead_i==1] <- 1
    
    d$wast_rec30d[d$wast_inc==1 & d$period_30d & d$rec_inc_lead_i==0 & d$wast_rec30d!=1] <- 0
    d$wast_rec60d[d$wast_inc==1 & d$period_60d & d$rec_inc_lead_i==0 & d$wast_rec60d!=1] <- 0
    d$wast_rec90d[d$wast_inc==1 & d$period_90d & d$rec_inc_lead_i==0 & d$wast_rec90d!=1] <- 0
    
    d$sevwast_inc30d[d$wast_inc==1 & d$period_30d & d$sev_inc_lead_i==0 & d$sevwast_inc30d!=1] <- 0
    d$sevwast_inc60d[d$wast_inc==1 & d$period_60d & d$sev_inc_lead_i==0 & d$sevwast_inc60d!=1] <- 0
    d$sevwast_inc90d[d$wast_inc==1 & d$period_90d & d$sev_inc_lead_i==0 & d$sevwast_inc90d!=1] <- 0
    
    #Stop for loop if all current leading observations are beyond 90 days
    if(min(d$days_lead_i, na.rm=T) > 90) break
  }
  
  #Drop intermediate variables
  d <- subset(d, select = -c(agelag, wastlag, sevwastlag, midpoint_age, wastchange, sevwastchange, past_wast, past_sevwast,
                             future_wast, future_sevwast,  sevwast_falter, sevwasting_episode_lag, sev_incident_age, sev_maxage,
                             sevduration, rec_inc_lead_i, sev_inc_lead_i, days_lead_i, period_30d,period_60d,period_90d, next_midpoint
  )) %>%
    ungroup() %>% as.data.frame()
  if(dropBornWasted==T){
    d <- subset(d, select = -c(born_wast_inc, born_sevwast_inc)) %>%
      ungroup() %>% as.data.frame()      
  }
  
  #merge back in other columns
  d <- merge(d, othercolumns, by=c("SUBJID", "AGEDAYS"))
  
  #rename columns to match other functions
  d <- d %>% rename(wast_rec = wast_rec_inc,
                    sevwast_rec = sevwast_rec_inc) 
  
  return(d)
}  





#----------------------------------------------
#Function to calculate summary tables
#----------------------------------------------

WastIncTable<-function(d, strat=T, agecats=c(6*30, 12*30, 18*30, 24*30), agecat_rownames=c("0-6 months","6-12 months", "12-18 months", "18-24 months")){
  
  
  if(strat==T){
    d$agecat <- as.factor(findInterval(d$AGEDAYS, agecats, rightmost.closed=F))
  }  
  
  summary<-WastIncSummary(d, strat=strat)
  tab<-summary[[1]]
  means<-data.frame(strata=rep("Overall",nrow(summary[[2]])), summary[[2]])
  
  if(strat==T){
    
    strattab<-stratmeans<-NULL
    for(i in 1:(length(agecats))){
      temp<-WastIncSummary(d[d$agecat==(i-1),], strat=F)
      strattab<-rbind(strattab, temp[[1]])
      stratmeans<-rbind(stratmeans, data.frame(strata=rep(agecat_rownames[i],nrow(temp[[2]])), temp[[2]]))
    }
    stratmeans<-stratmeans[!is.na(stratmeans$Mean),]
    
    #Replace truncated durations with stratified durations calculated from full data
    numcats<-length(unique(stratmeans$strata))
    #if(numcats>nrow(summary[[3]])){numcats<-nrow(summary[[3]])}
    strattab$total_duration[!is.na(strattab$total_duration)] <- summary[[3]]$total_duration[1:numcats]
    strattab$average_duration[!is.na(strattab$total_duration)]  <- summary[[3]]$average_duration[1:numcats]
    agelvls=levels(stratmeans$strat)
    for(i in 1:(numcats)){
      stratmeans[stratmeans$strata== agelvls[i] & 
                   stratmeans$statistic=="Average\nduration\nof\nwasting" ,3:8] <- summary[[4]][summary[[4]]$agecat==(i-1), -1]
    }
    
    
    tab<-rbind(tab, strattab)
    means<-rbind(means, stratmeans)
    rownames(tab)<-c("Overall",agecat_rownames)
    
  }
  
  tf<-tab_format(tab)
  
  
  
  return(list(tab1=tf[[1]], tab2=tf[[2]], tab3=tf[[3]], tab=tf[[4]], means=means))
  
  
  
}



#----------------------------------------------
#Function to calculate summary statistics
#----------------------------------------------

WastIncSummary<-function(d, strat=F){
  require(epitools)
  
  #Average number of measurements per child
  child_nmeas <- d %>% group_by(SUBJID) %>%
    summarize(num_measurements=n()) %>% ungroup() %>% 
    summarize(num_children=n(),
              num_measurements=mean(num_measurements))
  
  
  #Overall sum and means
  Incidence_df <- d %>% 
    summarize(      
      prev_wast= mean(WHZ < -2, na.rm=T) * 100,
      prev_sevwast= mean(WHZ < -3, na.rm=T) * 100,
      anywast= ifelse(sum(wast) > 0 ,1,0) * 100,
      anysevwast= ifelse(sum(sevwast) > 0 ,1,0) * 100,
      persontime=sum(wast_risk*period_length, na.rm=T),
      sev_persontime=sum(sevwast_risk*period_length, na.rm=T),
      recovery_persontime=sum(wast_rec_risk*period_length, na.rm=T),
      sevrecovery_persontime=sum(sevwast_risk*period_length, na.rm=T),
      wast_ep=sum(wast_inc, na.rm=T),
      sevwast_ep=sum(sevwast_inc, na.rm=T),
      wast_rec_ep=sum(wast_rec, na.rm=T),
      sevwast_rec_ep=sum(sevwast_rec, na.rm=T),
      recoveries30d=sum(wast_rec30d==1, na.rm=T),
      recoveries60d=sum(wast_rec60d==1, na.rm=T),
      recoveries90d=sum(wast_rec90d==1, na.rm=T),
      sevwast30d=sum(sevwast_inc30d==1, na.rm=T),
      sevwast60d=sum(sevwast_inc60d==1, na.rm=T),
      sevwast90d=sum(sevwast_inc90d==1, na.rm=T),
      no_recoveries30d=sum(wast_rec30d==0, na.rm=T),
      no_recoveries60d=sum(wast_rec60d==0, na.rm=T),
      no_recoveries90d=sum(wast_rec90d==0, na.rm=T),
      no_sevwast30d=sum(sevwast_inc30d==0, na.rm=T),
      no_sevwast60d=sum(sevwast_inc60d==0, na.rm=T),
      no_sevwast90d=sum(sevwast_inc90d==0, na.rm=T)
    ) %>% 
    mutate(  wastIR=wast_ep/persontime * 1000,
             sevwastIR=sevwast_ep/sev_persontime * 1000,
             wastrecIR=wast_rec_ep/recovery_persontime * 1000,
             sevwastrecIR=sevwast_rec_ep/sevrecovery_persontime * 1000,
             perc_wastrec_30d= sum(recoveries30d)/sum(no_recoveries30d)*100,
             perc_wastrec_60d= sum(recoveries60d)/sum(no_recoveries60d)*100,
             perc_wastrec_90d= sum(recoveries90d)/sum(no_recoveries90d)*100,
             perc_sevwastinc_30d= sum(sevwast30d)/sum(no_sevwast30d)*100,
             perc_sevwastinc_60d= sum(sevwast60d)/sum(no_sevwast60d)*100,
             perc_sevwastinc_90d= sum(sevwast90d)/sum(no_sevwast90d)*100
    ) %>%
    as.data.frame()
  
  #Calculate average episode lengths
  d <- d %>% group_by(SUBJID) %>%
    mutate(state_run = cumsum( wast_inc+sevwast_inc+wast_rec) + 1,
           wast_run = ifelse(wast==1, state_run, 0)) %>%
    ungroup() %>%  group_by(SUBJID, state_run) %>% 
    mutate(state_dur = sum(period_length),
           wast_dur = ifelse(wast==1, state_dur, 0)) %>%
    as.data.frame()
  
  #data frame of episode durations
  episode_duration <- d %>% 
    filter(wast_inc==1) %>% #drop non-wasting periods
    group_by(SUBJID, state_run ) %>% 
    slice(1) %>% ungroup() %>% 
    subset(., select=c(SUBJID,agecat,wast_dur)) %>% 
    as.data.frame()
  
  
  # #Calculate mean, max, and total duration of wasting per child
  duration <- episode_duration %>%
    group_by(SUBJID) %>%
    summarize(
      total_duration= sum(wast_dur, na.rm=T)
    )
  #Average episode length
  average_duration=mean(episode_duration$wast_dur)
  #average total time wasted for each child
  total_duration=mean(duration$total_duration)
  
  if(strat==T){
    duration_strat_average_duration <- episode_duration %>% group_by(agecat) %>% 
      summarize(average_duration=mean(wast_dur))
    duration_strat_total_duration <- episode_duration %>% group_by(agecat, SUBJID) %>% 
      summarize(total_duration=sum(wast_dur)) %>% 
      ungroup() %>% group_by(agecat) %>% 
      summarize(total_duration=mean(total_duration))
    duration_strat <- merge(duration_strat_average_duration, duration_strat_total_duration, by="agecat")
  }
  
  
  #Create taable of summary statistics
  tab <- data.frame(child_nmeas, Incidence_df[1:4], total_duration, average_duration, Incidence_df[-c(1:4, 13:24)])
  
  #Calculate means
  means <- NULL
  try(
    means <- rbind(
      #longitudinal prevalence
      mean95CI(Y=(d$WHZ < -2), id=d$SUBJID, proportion=T, percent=F),
      mean95CI(Y=(d$WHZ < -3), id=d$SUBJID, proportion=T, percent=F),
      #duration
      mean95CI(Y=episode_duration$wast_dur, id=episode_duration$SUBJID, proportion=F, percent=F),
      #incidence rates
      mean95CI(tab$wast_ep, persontime=tab$persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$sevwast_ep, persontime=tab$sev_persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$wast_rec_ep, persontime=tab$recovery_persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$sevwast_rec_ep, persontime=tab$sevrecovery_persontime, count=T) * c(1,rep(1000,5)),
      #percent recovery
      mean95CI(Y=d$wast_rec30,  proportion=T, percent=F),
      mean95CI(Y=d$wast_rec60,  proportion=T, percent=F),
      mean95CI(Y=d$wast_rec90,  proportion=T, percent=F),
      mean95CI(Y=d$sevwast_inc30d,  proportion=T, percent=F),
      mean95CI(Y=d$sevwast_inc60d,  proportion=T, percent=F),
      mean95CI(Y=d$sevwast_inc90d,  proportion=T, percent=F)))
  
  means <- data.frame(statistic =   c("Prevalence\nof\nwasting",
                                      "Prevalence\nof\nsevere\nwasting",
                                      "Average\nduration\nof\nwasting",
                                      "Wasting\nincidence\nrate",
                                      "Severe\nwasting\nincidence\nrate",
                                      "Wasting\nrecovery\nincidence\nrate",
                                      "Severe\nwasting\nrecovery\nincidence\nrate",
                                      "Percent\nwasting\nrecovered\nin 30 days",
                                      "Percent\nwasting\nrecovered\nin 60 days",
                                      "Percent\nwasting\nrecovered\nin 90 days",
                                      "Percent\nfalter to\nsevere\nwasting\nin 30 days",
                                      "Percent\nfalter to\nsevere\nwasting\nin 60 days",
                                      "Percent\nfalter to\nsevere\nwasting\nin 90 days"),means)
  
  if(strat==T){
    duration_strat_mean95CI <- NULL
    for(i in levels(d$agecat)){
      temp<-data.frame(agecat=i ,mean95CI(Y=episode_duration$wast_dur[episode_duration$agecat==i], id=episode_duration$SUBJID[episode_duration$agecat==i]))
      duration_strat_mean95CI <- rbind(duration_strat_mean95CI, temp)
    }
    
    return(list(tab, means, duration_strat, duration_strat_mean95CI))
  }else{
    return(list(tab, means))
  }
}









#----------------------------------------------
#Function to format summary tables
#----------------------------------------------
tab_format <- function(tab){
  
  
  tab<-as.data.frame(tab)
  tab$`Child age stratification` <- rownames(tab)
  tab1 <- tab %>% 
    subset(., select=c(
      `Child age stratification`,
      num_children,
      num_measurements,
      prev_wast,
      prev_sevwast,
      anywast,
      anysevwast,
      #total_duration,
      average_duration
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Prevalence of wasting across all measurements` =  prev_wast,
            `Prevalence of severe wasting across all measurements` = prev_sevwast,
            `Proportion of children who were ever wasted`= anywast,
            `Proportion of children who were ever severely wasted` = anysevwast,
            #`Total duration of wasting episodes (days)` = total_duration,
            `Average duration of wasting episodes (days)` = average_duration)
  
  tab2 <- tab %>% 
    subset(., select=c(
      `Child age stratification`,
      num_children,
      num_measurements,
      wast_ep,
      sevwast_ep,
      persontime,
      sev_persontime,
      wastIR,
      sevwastIR
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Number of wasting episodes` = wast_ep,
            `Number of severe wasting episodes` = sevwast_ep,
            `No. of days at risk of wasting` = persontime,
            `No. of days at risk of severe wasting` = sev_persontime,
            `Wasting incidence rate per 1000 days` = wastIR,
            `Severe wasting incidence rate per 1000 days` = sevwastIR
    )
  
  
  tab3 <- tab %>% 
    subset(., select=c(
      `Child age stratification`,
      num_children,
      num_measurements,
      wast_rec_ep,
      recovery_persontime,
      wastrecIR,
      perc_wastrec_30d,       
      perc_wastrec_60d,
      perc_wastrec_90d
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Number of recoveries from wasting` = wast_rec_ep,
            `No. of days at eligible for wasting recovery` = recovery_persontime,
            `Wasting recovery incidence rate per 1000 days` = wastrecIR,
            `Percent of wasting episodes recovered from in 30 days` = perc_wastrec_30d,       
            `Percent of wasting episodes recovered from in 60 days` = perc_wastrec_60d,
            `Percent of wasting episodes recovered from in 90 days` = perc_wastrec_90d)
  return(list(tab1,tab2,tab3, tab))
}






#----------------------------------------------
#Plot functions
#----------------------------------------------

#GAM curve
WHZ_curve<-function(df){
  theme_set(theme_bw())
  
  # grab a color blind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(1,3,7)]
  
  p <- gmsn_age<-ggplot(df, aes(x = AGEDAYS)) +
    geom_smooth(aes(y=WHZ), color="#D55E00") +
    geom_jitter(aes(y=WHZ,x=AGEDAYS), height = 0.2, width=0.2,  alpha = 0.1, size=0.5)+
    labs(y = "WHZ",
         x = "Child Age (Days)",
         title = "") +
    theme(strip.background = element_blank())
  
  return(p)
}


#spaghetti plot
spaghetti<-function(df){
  theme_set(theme_bw())
  
  set.seed(12345)
  d <- df %>% group_by(SUBJID) %>%
    mutate(alpha=  ifelse(runif(n = 1) > 0.95 , 1, 0.1))
  
  p <- ggplot(d, aes(x=AGEDAYS, y=WHZ)) + 
    geom_line() + guides(colour=FALSE) + xlab("Child Age (Days)") +
    ylab("WHZ") + aes(alpha=alpha, group=factor(SUBJID)) + guides(alpha=FALSE)
  
  return(p)
}


#Heatmap
heatmap<-function(d){
  theme_set(theme_bw())
  
  #generate whz categories
  d$wastcat <- 0
  d$wastcat[d$WHZ<(-2)] <- -1
  d$wastcat[d$WHZ<(-3)] <- -2
  
  #Sum an ad-hoc total wasting score to rank children by  
  d <- d %>% group_by(SUBJID) %>% mutate(wastscore=sum(wastcat, na.rm=T)) %>% ungroup()
  table(d$wastscore)
  
  #Add level for missingness
  d$wastcat[is.na(d$WHZ)] <- 1
  d$wastcat <- as.factor(d$wastcat)
  d$SUBJID <- as.factor(d$SUBJID)
  table(d$wastcat)
  
  #make ordered childnum by amount of wasting
  d <- d %>%
    arrange(wastscore) %>%  
    mutate(SUBJID = factor(SUBJID, unique(SUBJID))) 
  d$childnum <- as.numeric(d$SUBJID)
  head(as.data.frame(d))
  
  #Create a child age in months variable
  d$agemonths<-floor(d$AGEDAYS/30)
  
  
  if(sum(d$wastcat== -2) > 0){
    levels(d$wastcat)<-c("Severely wasted", "Wasted", "Not wasted")
    cbPalette <- c("#56B4E9", "#E69F00", "#c1c1c1")
  }else{
    levels(d$wastcat)<-c( "Wasted", "Not wasted")
    cbPalette <- c( "#E69F00", "#c1c1c1")   
  }
  
  p <- ggplot(d,aes(x=agemonths,y=childnum)) + 
    geom_tile(aes(fill=wastcat)) + scale_fill_manual(values=cbPalette) +
    xlab("Child Age (Months)") + ylab("Child number")
  ggtitle("Heatmap of wasting and severe wasting episodes") +
    theme(legend.title=element_blank())
  
  return(p)
}





#Incidence/prevalence figures
means_plot<-function(df){
  
  #df<-df[df$strata!=">24 months",]
  df$statistic <- factor(df$statistic, levels = unique(df$statistic))
  
  theme_set(theme_bw())
  cbPalette <- c( "#56B4E9" , "#999999" , "#999999", "#999999" , "#999999", "#999999" ) #, #f7a809, "#56B4E9",  "#E69F00",)
  p <- ggplot(df, aes(`Child age stratification`)) + 
    geom_point(aes(x=strata, y=Mean, fill=strata, color=strata), size = 4) +
    geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=strata), 
                   alpha=0.5, size = 3) +
    facet_wrap(~statistic, scales = "free")  + 
    scale_fill_manual(values=cbPalette) +
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.background = element_blank(),
          legend.position="none")
  
  return(p)
}










#Pooled logistic functions:

# --------------------------------------
# Robust clustered SE function
# http://people.su.se/~ma/mcluster.R
# R (www.r-project.org) codes for computing multi-way clustered-standard errors
# Mahmood Arai, Jan 21, 2008. 
# See: Thompson (2006), Cameron, Gelbach and Miller (2006) and Petersen (2006).
#
# slightly modified to have it return the vcovCL object
# rather than the updated fit (since might need the VC matrix)
# --------------------------------------
cl   <- function(dat,fm, cluster){
  # dat: data used to fit the model
  # fm : model fit (object)
  # cluster : vector of cluster IDs
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  return(vcovCL)
}

# --------------------------------------
# function to estimate poisson models
# with robust SEs
# --------------------------------------
poissonRB <- function(fmla,dat,print=TRUE) {
  # poisson regression with robust SEs
  # fmla : formula for the model fit
  # dat  : data used to fit the model (has to include "SUBJID" for individuals)
  # print: print results ?
  
  # restrict to complete cases
  dat <- dat[complete.cases(dat[,c(all.vars(paste0(fmla)))]),] 
  
  fit <- glm(fmla,family=poisson,data=dat,model=FALSE,x=TRUE)
  
  dat <- na.omit(dat[ , c("SUBJID", all.vars(formula(fit)))])
  vcovCL <- cl(dat=dat,fm=fit,cluster=dat$SUBJID)
  rfit <- coeftest(fit, vcovCL)
  if(print==TRUE) {
    cat(paste("N obs=",nrow(dat)))
    print(rfit)
  }
  list(fit=fit,rfit=rfit,vcovCL=vcovCL)
}


# --------------------------------------
# Wrapper function to calculate RR of 
# any wasting by 6 months and recovery in 60
# days among those wasted using TMLE, and "HR"/intensity
# ratio of incident wasting using pooled logistic regressions
# --------------------------------------

riskfactor4b<-function(df, 
                       A="BIRTHWT",
                       Wvars=NULL, 
                       n.cat=4, 
                       reflevel=1, 
                       SLlibrary="SL.glm", 
                       Acuts=NULL, 
                       Alevels=NULL, 
                       agerange=c(0,6),
                       born_not_wast=F,
                       overall_dist=T,
                       run_anywast=T,
                       run_recovery=T,
                       run_irr=T){
  
  cat("\n", unique(df$STUDYID),": ", unique(df$COUNTRY),"\n")
  
  Avar<-A
  
  #Rename risk factor to A
  colnames(df)[grep(paste0("^",A,"$"), colnames(df))]<-"A"
  
  #Drop any missing A
  df <- df %>% filter(!is.na(A))
  
  if(born_not_wast==T){
    birth_df<-df %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
      slice(which.min(AGEDAYS)) %>% 
      mutate(bornwast=as.numeric(WHZ < -2))
    birth_df$bornwast[birth_df$AGEDAYS > 60] <- NA
    
    df <- left_join(df, birth_df, by= c("STUDYID", "SUBJID"))
    df <- df %>% filter(bornwast==0)
  }
  
  #Tabulate n in agerange (to diagnose when no estimates are returned)
  n_in_agerange<- nrow(df[df$AGEDAYS>=agerange[1] & df$AGEDAYS<agerange[2],])
  
  inc <- df %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
    mutate(wastcount=cumsum(wast_inc),
           anywast=ifelse(wastcount>0,1,0)) %>%
    slice(which.min(abs(AGEDAYS - agerange[2]*30))) %>% #Find observation with age closest to selected age
    filter((AGEDAYS > (agerange[2]-2)*30) ) %>% # Last observation must be within 2 months of the age of cumulative incidence 
    ungroup() %>%
    as.data.frame()
  
  #Calculate the average and min age used in the cumulative inc. calculations
  if(nrow(inc)>0){
    mean_age<-mean(inc$AGEDAYS)
    min_age<-min(inc$AGEDAYS)
  }else{
    mean_age<-0
    min_age<-0
  }
  
  if(is.null(Wvars)){
    inc$W1<-rep(1, nrow(inc))
    inc$W2<-rep(1, nrow(inc))
    Wvars<-c("W1","W2")
  }
  
  
  #Set acuts based on individual distribution if needed
  if(is.null(Acuts) | overall_dist==F){
    Acuts=quantile(inc$A, probs = c((1:(n.cat-1))/n.cat), na.rm=T)
  }
  
  if(run_anywast==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),": ",dim(df),"\n")
    Anywast6mo<-NULL
    try(
      Anywast6mo<-tmle_risk(dat=as.data.frame(inc),
                            W=Wvars,
                            A="A",
                            n.cat=n.cat,
                            reflevel=reflevel,
                            Acuts=Acuts,
                            Alevels=Alevels,
                            outputdf=NULL,
                            Y="anywast",
                            family="binomial",
                            SLlibrary=SLlibrary,
                            overall.dist=overall_dist,
                            sparseN=0,
                            adjusted=F))
    try(
      Anywast6mo<-data.frame(study=rep(unique(df$STUDYID), nrow(Anywast6mo)) ,country=rep(unique(df$COUNTRY), nrow(Anywast6mo)) , variable=rep(Avar,  nrow(Anywast6mo)),
                             Anywast6mo, n_in_agerange=rep(n_in_agerange, nrow(Anywast6mo))
                             , mean_age=rep(mean_age, nrow(Anywast6mo))
                             , min_age=rep(min_age, nrow(Anywast6mo))))
    
    if(is.null(Anywast6mo)){
      Anywast6mo<-data.frame(study=unique(df$STUDYID),country=unique(df$COUNTRY), variable=Avar,
                             variable.1=NA, level=NA,         ATE=NA,
                             ATE.var=NA,     ATE.CI1=NA,      ATE.CI2=NA,     ATE.Pval=NA,        RR=NA,    RR.CI1=NA,
                             RR.CI2=NA,      RR.Pval=NA,   logRR.psi=NA,   logRR.var=NA, compN=NA, refN=NA,    a=NA,    b=NA,
                             c=NA,    d=NA, meanLevel=NA, meanN=NA,     meanY=NA,   mean.sd=NA,     
                             mean.se=NA,  mean.CI1=NA, mean.CI2=NA, n_in_agerange=NA, mean_age=mean_age, min_age=min_age)
    }
  }else{
    Anywast6mo<-NULL
  }
  
  
  if(run_recovery==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),"\n")
    drec<- df %>% filter(wast_inc==1 & !is.na(A)  & AGEDAYS >= agerange[1]*30 & AGEDAYS < agerange[2]*30) %>%
      as.data.frame()
    
    if(Wvars[1]=="W1"){
      drec$W1<-rep(1, nrow(drec))
      drec$W2<-rep(1, nrow(drec))
    }
    
    #Need to add in SUBJID as clustering variable
    wastrec60d<-NULL
    try(wastrec60d<-tmle_risk(dat=as.data.frame(drec),
                              W=Wvars,
                              A="A",
                              n.cat=n.cat,
                              reflevel=reflevel,
                              Acuts=Acuts,
                              Alevels=Alevels,
                              outputdf=NULL,
                              Y="wast_rec60d",
                              family="binomial",
                              SLlibrary="SL.glm",
                              overall.dist=overall_dist,
                              sparseN=0,
                              adjusted=F))
    try(
      wastrec60d<-data.frame(study=rep(unique(df$STUDYID), nrow(wastrec60d)) , country=rep(unique(df$COUNTRY), nrow(wastrec60d)) , variable=rep(Avar,  nrow(wastrec60d)),
                             wastrec60d, n_in_agerange=rep(n_in_agerange, nrow(wastrec60d))))
    
    if(is.null(wastrec60d)){
      wastrec60d<-data.frame(study=unique(df$STUDYID),country=unique(df$COUNTRY), variable=Avar,
                             variable.1=NA, level=NA,         ATE=NA,
                             ATE.var=NA,     ATE.CI1=NA,      ATE.CI2=NA,     ATE.Pval=NA,        RR=NA,    RR.CI1=NA,
                             RR.CI2=NA,      RR.Pval=NA,   logRR.psi=NA,   logRR.var=NA, compN=NA, refN=NA,    a=NA,    b=NA,
                             c=NA,    d=NA, meanLevel=NA, meanN=NA,     meanY=NA,   mean.sd=NA,     
                             mean.se=NA,  mean.CI1=NA, mean.CI2=NA, n_in_agerange=NA)
    }
    
  }else{
    wastrec60d<-NULL
  }
  
  
  
  #pooled regression for wasting intensity
  if(run_irr==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),"\n")
    wastrisk <- df %>% filter(wast_risk==1 & !is.na(SUBJID)  & !is.na(AGEDAYS)  & !is.na(wast)  & AGEDAYS >= agerange[1]*30 & AGEDAYS < agerange[2]*30) %>% arrange(A)
    if(overall_dist==F){
      Acuts<-quantile(wastrisk$A, probs = c((1:(n.cat-1))/n.cat), na.rm=T)
    }
    wastrisk$Acat <- factor(findInterval(wastrisk$A, Acuts, left.open=T))
    wastrisk$Acat <- relevel(wastrisk$Acat, ref=levels(wastrisk$Acat)[reflevel])
    wastrisk$ageZ<- scale(wastrisk$AGEDAYS, center = TRUE, scale = TRUE)
    irr <- poissonRB(fmla="wast~Acat + poly(ageZ,3,raw=TRUE)", dat=wastrisk)
    logres<-irr$rfit
    
    res<-data.frame(study=rep(unique(df$STUDYID), nrow(logres)),country=rep(unique(df$COUNTRY), nrow(logres)) , variable=rep(Avar,  nrow(logres)), coef=rownames(logres), RR=exp(logres[,1]), ci.lb= exp(logres[,1]-1.96*logres[,2]), ci.ub=  exp(logres[,1]+1.96*logres[,2]),logRR=logres[,1] , log.se=logres[,2], pvalue=logres[,4])
  }else{
    res<-NULL
  }
  return(list(Anywast6mo=Anywast6mo, wastrec60d=wastrec60d, irr=res))
}












# WastIncCalc_old<-function(d, washout=60){
#   require(tidyverse)
#   require(zoo)  
#   
#   
#   
#   #Filter out extreme or missing whz values
#   d <- d %>%  filter(!is.na(WHZ)) %>%
#     filter(WHZ > (-5) & WHZ < 5)
#   
#   #Remove duplicate ages
#   ndropped <- nrow(d[duplicated(cbind(d$SUBJID, d$AGEDAYS)), ])
#   d <- d[!duplicated(cbind(d$SUBJID, d$AGEDAYS)), ]
#   if(ndropped>0) cat("\n",ndropped," observations dropped due to duplicate ages\n")
#   
#   
#   #Create visit variable
#   #(This replaces the visit variables, but some studies have missing 
#   #visit data)
#     d <- d %>% group_by(SUBJID) %>% mutate(VISITNUM = rank(AGEDAYS))
#   
#   
#   #Extract required columns and save others to merge back in later
#   othercolumns <- d %>% subset(., select= -c(WHZ, VISITNUM)) 
#   d <- d %>% subset(., select= c(SUBJID, WHZ, AGEDAYS, VISITNUM)) 
#   
#   
#   #generate wasting
#   d$wast= ifelse(d$WHZ < (-2),1,0)
#   d$sevwast= ifelse(d$WHZ < (-3),1,0)
#   
#   d <- d %>%
#     arrange(SUBJID, AGEDAYS) %>%
#     group_by(SUBJID) %>%
#     mutate(wastlag=lag(wast),
#            sevwastlag=lag(sevwast),
#            agelag=lag(AGEDAYS),
#            agelead=lead(AGEDAYS),
#            period_lag= AGEDAYS - agelag, #Calculate number of days since last measurement
#            period_lead= agelead - AGEDAYS,
#            wastchange= as.numeric((wast - wastlag)!=0 & !is.na(wastlag)),
#            sevwastchange=  as.numeric(((sevwast - sevwastlag)==1 & !is.na(sevwastlag)))) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   #Fill in missing lag variables with 0 (for the first measurements)
#   d$wastlag[is.na(d$wastlag)] <- 0
#   d$sevwastlag[is.na(d$sevwastlag)] <- 0
#   d$agelag[is.na(d$agelag)] <- 0
#   d$period_lag[is.na(d$period_lag)] <- d$AGEDAYS[is.na(d$period_lag)] #Set days since birth
#   d$period_lead[is.na(d$period_lead)] <- 30  #Assume month if no further observations
#   d$period_length<-(d$period_lag + d$period_lead)/2 #Assumed period length
#   
#   #Expand longform dataset to have a row for every child day
#   vals <- expand.grid(AGEDAYS = 1:max(d$AGEDAYS),
#                       SUBJID = unique(d$SUBJID))
#   d <- left_join(vals, d, by=c("SUBJID", "AGEDAYS"))
#   dfill<-d %>% subset(.,select=c(AGEDAYS,SUBJID,wast, sevwast)) %>% 
#     group_by(SUBJID) %>%
#     na.locf(., fromLast=F) %>%
#     rename(wastfill=wast, sevwastfill=sevwast) %>% 
#     arrange(AGEDAYS) %>% ungroup()
#   d <- merge(d, dfill, by=c("SUBJID","AGEDAYS"))
#   d <- d %>% arrange(SUBJID,AGEDAYS)
#   
#   d <- d %>%
#     arrange(SUBJID,AGEDAYS) %>%
#     group_by(SUBJID) %>%
#     mutate(
#       dayswast2m = roll_sum_fun(wastfill, washout),
#       dayssevwast2m = roll_sum_fun(sevwastfill, washout),
#       pastdayswast2m = lag_sum_fun(wastfill, washout),
#       pastdayssevwast2m = lag_sum_fun(sevwastfill, washout)) %>%
#     filter(!is.na(WHZ)) %>%
#     mutate( #Set past wasting to 0 for first observation
#       pastdayswast2m = ifelse(VISITNUM==1, 0, pastdayswast2m),
#       pastdayssevwast2m = ifelse(VISITNUM==1, 0, pastdayssevwast2m),
#       future_wast = ifelse(dayswast2m>0, 1,0),
#       future_sevwast = ifelse(dayssevwast2m>0, 1,0),
#       past_wast = ifelse(pastdayswast2m>0, 1,0),
#       past_sevwast = ifelse(pastdayssevwast2m>0, 1,0)) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   #Calculate incidence
#   # Severe recovery incidence requires first full recovery from wasting
#   d <- d %>% group_by(SUBJID) %>%
#     mutate(
#       wast_risk= as.numeric(wastlag==0 & past_wast==0),
#       sevwast_risk= as.numeric(sevwastlag==0 & past_sevwast==0),
#       wast_rec_risk= as.numeric(past_wast==1),
#       sevwast_rec_risk= as.numeric(past_sevwast==1),
#       wast_inc= as.numeric(wast==1 & wastlag==0 & past_wast==0),
#       sevwast_inc= as.numeric(sevwast==1 & sevwastlag==0 & past_sevwast==0),
#       wast_rec= as.numeric(wast==0 & wastlag==1 & future_wast==0), 
#       times_wasted=cumsum(wast_inc),
#       times_sevwasted=cumsum(sevwast_inc),
#       times_rec=cumsum(wast_rec),
#       times_rec_lag=lag(times_rec),
#       sevwast_rec=ifelse(times_rec!=times_rec_lag & times_rec==times_sevwasted,1,0)
#     ) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   
#   #Calculate 30,60, 90 day recovery
#   #And 30, 60, 90 day faltering into severe wasting
#   #Expand longform dataset to have a row for every child day
#   dsub <- d %>% subset(.,select=c(SUBJID, AGEDAYS, wast_inc, wast_rec, sevwast_inc))
#   d_expand <- left_join(vals, dsub, by=c("SUBJID", "AGEDAYS"))
#   d_rec <- d_expand %>% group_by(SUBJID) %>% arrange(SUBJID,AGEDAYS) %>%
#     mutate(
#       wast_rec30d= as.numeric(roll_sum_fun(wast_rec, 30)>0 & (wast_inc==1 & !is.na(wast_inc))),
#       wast_rec60d= as.numeric(roll_sum_fun(wast_rec, 60)>0 & (wast_inc==1 & !is.na(wast_inc))),
#       wast_rec90d= as.numeric(roll_sum_fun(wast_rec, 90)>0 & (wast_inc==1 & !is.na(wast_inc))),     
#       
#       sevwast_inc30d= as.numeric(roll_sum_fun(sevwast_inc, 30)>0 & (wast_inc==1 & !is.na(wast_inc))),
#       sevwast_inc60d= as.numeric(roll_sum_fun(sevwast_inc, 60)>0 & (wast_inc==1 & !is.na(wast_inc))),
#       sevwast_inc90d= as.numeric(roll_sum_fun(sevwast_inc, 90)>0 & (wast_inc==1 & !is.na(wast_inc)))   
#     ) %>% 
#     filter(!is.na(wast_inc)) %>%
#     subset(.,select=c(SUBJID, AGEDAYS, wast_rec30d, wast_rec60d, wast_rec90d, sevwast_inc30d, sevwast_inc60d, sevwast_inc90d)) %>%
#     ungroup() %>% as.data.frame()
#   
#   
#   d <- merge(d, d_rec, by=c("SUBJID","AGEDAYS"))
#   d <- d %>% arrange(SUBJID,AGEDAYS)
#   
#   #Drop intermediate variables
#   d <- d %>% subset(., select= -c(wastlag, sevwastlag, agelag, agelead, 
#                                   period_lag, period_lead, wastchange, sevwastchange, 
#                                   wastfill, sevwastfill, dayswast2m, dayssevwast2m, pastdayswast2m,
#                                   pastdayssevwast2m, future_wast, future_sevwast, past_wast, past_sevwast))
#   #merge back in other columns
#   d <- merge(d, othercolumns, by=c("SUBJID", "AGEDAYS"))
#   
#   return(d)
# }
# 


# 
# #----------------------------------------------
# #Function to calculate summary statistics
# #----------------------------------------------
# 
# WastIncSummary_old<-function(d, strat=F){
#   require(epitools)
#   
#   
#   #Calculated days at risk per person
#   Incidence_df <- d %>% group_by(SUBJID) %>%
#     summarize(persontime=sum(wast_risk*period_length, na.rm=T),
#               sev_persontime=sum(sevwast_risk*period_length, na.rm=T),
#               recovery_persontime=sum(wast_rec_risk*period_length, na.rm=T),
#               sevrecovery_persontime=sum(sevwast_risk*period_length, na.rm=T),
#               wast_ep=sum(wast_inc, na.rm=T),
#               sevwast_ep=sum(sevwast_inc, na.rm=T),
#               wast_rec_ep=sum(wast_rec, na.rm=T),
#               sevwast_rec_ep=sum(sevwast_rec, na.rm=T),
#               recoveries30d=sum(wast_rec30d==1, na.rm=T),
#               recoveries60d=sum(wast_rec60d==1, na.rm=T),
#               recoveries90d=sum(wast_rec90d==1, na.rm=T),
#               sevwast30d=sum(sevwast_inc30d==1, na.rm=T),
#               sevwast60d=sum(sevwast_inc60d==1, na.rm=T),
#               sevwast90d=sum(sevwast_inc90d==1, na.rm=T),
#               no_recoveries30d=sum(wast_rec30d==0, na.rm=T),
#               no_recoveries60d=sum(wast_rec60d==0, na.rm=T),
#               no_recoveries90d=sum(wast_rec90d==0, na.rm=T),
#               no_sevwast30d=sum(sevwast_inc30d==0, na.rm=T),
#               no_sevwast60d=sum(sevwast_inc60d==0, na.rm=T),
#               no_sevwast90d=sum(sevwast_inc90d==0, na.rm=T)
#     ) %>% 
#     mutate(  wastIR=wast_ep/persontime * 1000,
#              sevwastIR=sevwast_ep/sev_persontime * 1000,
#              wastrecIR=wast_rec_ep/recovery_persontime * 1000,
#              sevwastrecIR=sevwast_rec_ep/sevrecovery_persontime * 1000,
#              perc_wastrec_30d= sum(recoveries30d)/sum(no_recoveries30d)*100,
#              perc_wastrec_60d= sum(recoveries60d)/sum(no_recoveries60d)*100,
#              perc_wastrec_90d= sum(recoveries90d)/sum(no_recoveries90d)*100,
#              perc_sevwastinc_30d= sum(sevwast30d)/sum(no_sevwast30d)*100,
#              perc_sevwastinc_60d= sum(sevwast60d)/sum(no_sevwast60d)*100,
#              perc_sevwastinc_90d= sum(sevwast90d)/sum(no_sevwast90d)*100
#     ) %>%
#     as.data.frame()
#   #Replace perc_wastrec with NA if no wasting occurred
#   Incidence_df$perc_wastrec_30d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df$perc_wastrec_60d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df$perc_wastrec_90d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df$perc_sevwastinc_30d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df$perc_sevwastinc_60d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df$perc_sevwastinc_90d[Incidence_df$wast_ep==0]<-NA
#   Incidence_df
#   
#   
#   
#   #Calculate average episode lengths
#   d <- d %>% group_by(SUBJID) %>%
#     mutate(state_run = cumsum( wast_inc+sevwast_inc+wast_rec) + 1,
#            wast_run = ifelse(wast==1, state_run, 0)) %>%
#     as.data.frame()
#   
#   d <- d %>% group_by(SUBJID, state_run) %>% 
#     mutate(state_dur = sum(period_length),
#            wast_dur = ifelse(wast==1, state_dur, 0)) %>%
#     as.data.frame()
#   
#   head(d,30)
#   
#   
#   #Calculate mean, max, and total duration of wasting per child
#   duration <- d %>% 
#     filter(wast_run!=0) %>% #drop non-wasting periods
#     group_by(SUBJID, state_run ) %>% 
#     slice(1) %>% #Grab first row (so no duplication of period days)
#     ungroup() %>%
#     group_by(SUBJID) %>%
#     summarize(
#       total_duration= sum(wast_dur, na.rm=T),
#       average_duration= mean(wast_dur, na.rm=T))
#   head(duration)
#   
#   if(strat==T){
#     duration_strat <- d %>% 
#       filter(wast_run!=0) %>% #drop non-wasting periods
#       group_by(SUBJID, state_run ) %>% 
#       slice(1) %>% #Grab first row (so no duplication of period days)
#       ungroup() %>%
#       group_by(SUBJID, agecat) %>%
#       summarize(
#         total_duration= sum(wast_dur, na.rm=T),
#         average_duration= mean(wast_dur, na.rm=T)) %>%
#       ungroup()
#     
#     duration_strat_mean <- duration_strat %>%
#       group_by(agecat) %>%
#       summarize(
#         total_duration= mean(total_duration, na.rm=T),
#         average_duration= mean(average_duration, na.rm=T))      
#   }
#   
#   
#   #summarize longitudinal and cumulative prevalence
#   sumdf<- d %>% group_by(SUBJID) %>%
#     summarize(
#       num_measurements=n(),
#       prev_wast= mean(WHZ < -2, na.rm=T),
#       prev_sevwast= mean(WHZ < -3, na.rm=T),
#       anywast= ifelse(sum(wast) > 0 ,1,0),
#       anysevwast= ifelse(sum(sevwast) > 0 ,1,0)
#     ) %>% as.data.frame()
#   head(sumdf)
#   
#   
#   
#   #merge in duration
#   sumdf <- left_join(sumdf, duration, by="SUBJID")
#   
#   #merge in incidence 
#   sumdf <- left_join(sumdf, Incidence_df, by="SUBJID")
#   
#   head(sumdf)
#   
#   #Calculate overall and age-stratified estimates
#   studymeans=as.data.frame(t(colMeans(sumdf[-1], na.rm =T)))
#   studysum=as.data.frame(t(colSums(sumdf[-1], na.rm =T)))
#   
#   #Replace IR calculations with summed episodes and persontime (rather than mean of individual IR)
#   studymeans$wastIR<-mean(sum(sumdf$wast_ep, na.rm=T)/sum(sumdf$persontime), na.rm=T)*1000
#   studymeans$sevwastIR<-mean(sum(sumdf$sevwast_ep, na.rm=T)/sum(sumdf$sev_persontime), na.rm=T)*1000
#   studymeans$wastrecIR<-mean(sum(sumdf$wast_rec_ep, na.rm=T)/sum(sumdf$recovery_persontime), na.rm=T)*1000
#   studymeans$sevwastrecIR<-mean(sum(sumdf$sevwast_rec_ep, na.rm=T)/sum(sumdf$sevrecovery_persontime), na.rm=T)*1000
#   
#   
# 
#   tab<-cbind(nrow(sumdf) ,studymeans[,c(1:7,22:31)], studysum[,c(8:21)])
#   names(tab)[1]<-"num_children"
#   tab[3:6]<-tab[3:6]*100
# 
#   wastIR.CI <- pois.exact(tab$wast_ep, pt = tab$persontime, conf.level = 0.95)[3:5] * 1000
#   sevwastIR.CI <- pois.exact(tab$sevwast_ep, pt = tab$sev_persontime, conf.level = 0.95)[3:5] * 1000
#   wastrecIR.CI <- pois.exact(tab$wast_rec_ep, pt = tab$recovery_persontime, conf.level = 0.95)[3:5] * 1000
#   sevwastrecIR.CI <- pois.exact(tab$sevwast_rec_ep, pt = tab$sevrecovery_persontime, conf.level = 0.95)[3:5] * 1000
#   
#   mean_wast_incIR <- data.frame(N=NA, Mean=wastIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=wastIR.CI[2] ,  Upper.95.CI=wastIR.CI[3] )
#   mean_sevwast_inc <- data.frame(N=NA, Mean=sevwastIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=sevwastIR.CI[2] ,  Upper.95.CI=sevwastIR.CI[3] )
#   mean_wastrecIR <- data.frame(N=NA, Mean=wastrecIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=wastrecIR.CI[2] ,  Upper.95.CI=wastrecIR.CI[3] )
#   mean_sevwastrecIR <- data.frame(N=NA, Mean=wastrecIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=wastrecIR.CI[2] ,  Upper.95.CI=wastrecIR.CI[3] )
#   
#   colnames(mean_wast_incIR) <- colnames(mean_sevwast_inc) <- colnames(mean_wastrecIR) <- colnames(mean_sevwastrecIR) <- c("N","Mean","SD","Robust SE", "Lower 95%CI", "Upper 95%CI") 
#   
#   #Calculate means
#   means <- NULL
#   try(
#     means <- rbind(
#       mean95CI(Y=sumdf$prev_wast, id=sumdf$SUBJID, proportion=T, percent=F),
#       mean95CI(Y=sumdf$prev_sevwast, id=sumdf$SUBJID, proportion=T, percent=F),
#       mean95CI(Y=sumdf$average_duration, id=sumdf$SUBJID),
#       
#       mean_wast_incIR,
#       mean_sevwast_inc,
#       mean_wastrecIR,
#       mean_sevwastrecIR,
#       
#       mean95CI(Y=sumdf$perc_wastrec_30d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_wastrec_60d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_wastrec_90d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevwastinc_30d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevwastinc_60d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevwastinc_90d, id=sumdf$SUBJID, proportion=T, percent=T)))
#   
#   means <- data.frame(statistic =   c("Prevalence\nof\nwasting",
#                                       "Prevalence\nof\nsevere\nwasting",
#                                       "Average\nduration\nof\nwasting",
#                                       "Wasting\nincidence\nrate",
#                                       "Severe\nwasting\nincidence\nrate",
#                                       "Wasting\nrecovery\nincidence\nrate",
#                                       "Severe\nwasting\nrecovery\nincidence\nrate",
#                                       "Percent\nwasting\nrecovered\nin 30 days",
#                                       "Percent\nwasting\nrecovered\nin 60 days",
#                                       "Percent\nwasting\nrecovered\nin 90 days",
#                                       "Percent\nfalter to\nsevere\nwasting\nin 30 days",
#                                       "Percent\nfalter to\nsevere\nwasting\nin 60 days",
#                                       "Percent\nfalter to\nsevere\nwasting\nin 90 days"),means)
#   
#   if(strat==T){
#     duration_strat_mean95CI <- NULL
#     for(i in levels(d$agecat)){
#       temp<-data.frame(agecat=i ,mean95CI(Y=duration_strat$average_duration[d$agecat==i], id=duration_strat$SUBJID[d$agecat==i]))
#       duration_strat_mean95CI <- rbind(duration_strat_mean95CI, temp)
#     }
#     
#     mean95CI(Y=sumdf$average_duration, id=sumdf$SUBJID)
#     
#     return(list(tab, means, duration_strat_mean, duration_strat_mean95CI))
#   }else{
#     return(list(tab, means))
#   }
# }