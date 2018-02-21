

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
# HAZ<-c(0,-3.5,0,-3.5,0,0,0,-3.5,0,
#        -2.5,-2.5,-2.5,-2.5,-2.5,
#        -2.5,0,-2.5,0,0)
# HAZ<-rnorm(3000,-2,1)
# d<-data.frame(SUBJID,AGEDAYS,HAZ)
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
# test<- stuntIncCalc(d,washout=60)
# test2<-stuntIncTable(test,
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
stuntIncCalc<-function(d, washout=60, dropBornstunted=F){
  require(tidyverse)
  require(zoo)  
  
  #Filter out extreme or missing HAZ values
  d <- d %>%  filter(!is.na(HAZ)) %>%
    filter(HAZ > (-5) & HAZ < 5)
  
  #Remove duplicate ages
  ndropped <- nrow(d[duplicated(cbind(d$SUBJID, d$AGEDAYS)), ])
  d <- d[!duplicated(cbind(d$SUBJID, d$AGEDAYS)), ]
  if(ndropped>0) cat("\n-----------------------------------\n",ndropped," observations dropped due to duplicate ages\n_----------------------------------\n")
  
  
  #Create visit variable
  #(This replaces the visit variables, but some studies have missing 
  #visit data)
  d <- d %>% group_by(SUBJID) %>% mutate(VISITNUM = rank(AGEDAYS))
  
  
  #Extract required columns and save others to merge back in later
  othercolumns <- d %>% subset(., select= -c(HAZ, VISITNUM)) 
  d <- d %>% subset(., select= c(SUBJID, HAZ, AGEDAYS, VISITNUM)) 
  
  
  #generate stunting and severe stunting indicators
  d$stunt= ifelse(d$HAZ < (-2),1,0)
  d$sevstunt= ifelse(d$HAZ < (-3),1,0)
  
  # #Generate variables for length of period in days between prior observation and current observations
  # #and the next observations and current observations. Also generate variables for if child changed from
  # #not stunted to stunted (or severe stunted) between the last observation and the current observation.
  d <- d %>%
    arrange(SUBJID, AGEDAYS) %>%
    group_by(SUBJID) %>%
    mutate(
      agelag=lag(AGEDAYS),
      stuntlag=lag(stunt),
      sevstuntlag=lag(sevstunt),
      midpoint_age = AGEDAYS - (AGEDAYS - agelag)/2,
      stuntchange = stunt - lag(stunt),
      sevstuntchange = sevstunt - lag(sevstunt)
    ) %>%
    as.data.frame()
  d$agelag[is.na(d$agelag)] <- 0
  d$stuntlag[is.na(d$stuntlag)] <- 0
  d$sevstuntlag[is.na(d$sevstuntlag)] <- 0
  d$stuntchange[is.na(d$stuntchange)] <- d$stunt[is.na(d$stuntchange)]
  d$sevstuntchange[is.na(d$sevstuntchange)] <- d$sevstunt[is.na(d$sevstuntchange)]
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
  d$future_sevstunt <- d$future_stunt <-  d$past_sevstunt <- d$past_stunt <- rep(0, N)
  
  for(i in 1:washout){
    
    d <- d %>% group_by(SUBJID) %>%
      mutate(
        stunt_lag_i  = lag(stunt, i),
        sevstunt_lag_i  = lag(sevstunt, i),
        days_lag_i = abs(lag(midpoint_age, i) - midpoint_age),
        stunt_lead_i  = lead(stunt, i),
        sevstunt_lead_i  = lead(sevstunt, i),
        days_lead_i = abs(lead(midpoint_age, i) - midpoint_age)
      )
    
    d$washout_period_lag[d$days_lag_i > washout] <- F
    d$washout_period_lead[d$days_lead_i > washout] <- F
    
    d$past_stunt[d$stunt_lag_i==1 & d$washout_period_lag==T] <- 1
    d$past_sevstunt[d$sevstunt_lag_i==1 & d$washout_period_lag==T] <- 1              
    d$future_stunt[d$stunt_lead_i==1 & d$washout_period_lead==T] <- 1
    d$future_sevstunt[d$sevstunt_lead_i==1 & d$washout_period_lead==T] <- 1
    
    #Stop for loop if all current leading and lagging observations are beyond washout period
    if(min(d$days_lead_i, na.rm=T) & min(d$days_lag_i, na.rm=T) > washout) break
  }
  
  d <- d %>% 
    subset(., select= -c(washout_period_lag, washout_period_lead, 
                         stunt_lag_i, sevstunt_lag_i, days_lag_i, 
                         stunt_lead_i, sevstunt_lead_i, days_lead_i)) %>% 
    ungroup() %>% as.data.frame()
  
  head(d,30)
  
  #---------------------------------------------------------
  #Calculate stunting and stunting recovery incidence and risk
  #---------------------------------------------------------
  d$stunt_rec_inc <- d$sevstunt_rec_inc <- d$stunt_inc <- d$sevstunt_inc <- rep(0, N)
  d$sevstunt_rec_risk <- d$stunt_rec_risk <- d$sevstunt_risk <- d$stunt_risk <-  rep(0, N)
  
  
  d$stunt_inc[d$stuntchange==1 & d$past_stunt==0] <- 1 #stunting incidence if at risk of stunting and change in status between prior and current observation
  d$stunt_rec_inc[d$stuntchange== -1 & d$future_stunt==0] <- 1 #Recovery from stunting if status change to not stunted and no new stunting in the future washout period
  
  d$stunt_inc[d$stuntchange!=1 | d$past_stunt==1] <- 0 
  d$stunt_rec_inc[d$stuntchange!= -1 | d$future_stunt==1] <- 0  
  
  #Remove incidences of stunting if there has not been recovery from prior stunting episode
  #Is there a cleaner way of preventing recording the incidences earlier?
  table(d$stunt_inc)
  
  d <- d %>% group_by(SUBJID) %>% 
    mutate(sum_stunt_inc=cumsum(stunt_inc),
           sum_stunt_rec=cumsum(stunt_rec_inc))
  for(i in 1:nrow(d)){
    if(d$stunt_inc[i]==1 & (d$sum_stunt_inc[i]-d$sum_stunt_rec[i] > 1)){
      d$stunt_inc[i] <- 0 
      d <- d %>% group_by(SUBJID) %>% 
        mutate(sum_stunt_inc=cumsum(stunt_inc),
               sum_stunt_rec=cumsum(stunt_rec_inc))        
    }
    if(d$stunt_rec_inc[i]==1 & (d$sum_stunt_inc[i]-d$sum_stunt_rec[i] < 0)){
      d$stunt_rec_inc[i] <- 0 
      d <- d %>% group_by(SUBJID) %>% 
        mutate(sum_stunt_inc=cumsum(stunt_inc),
               sum_stunt_rec=cumsum(stunt_rec_inc))  
    }
  }
  table(d$stunt_inc)
  
  
  d <- subset(d, select = -c(sum_stunt_inc,sum_stunt_rec))
  #Make sure there isn't double recovery
  # d <- d %>% group_by(SUBJID) %>% 
  #   mutate(sum_stunt_inc=cumsum(stunt_inc),
  #          sum_stunt_rec=cumsum(stunt_rec_inc))
  # table(d$sum_stunt_inc- d$sum_stunt_rec != 1 & d$sum_stunt_inc- d$sum_stunt_rec != 0)
  # d$flag <- d$sum_stunt_inc- d$sum_stunt_rec != 1 & d$sum_stunt_inc- d$sum_stunt_rec != 0
  # 
  #Indicate length of incident episodes
  d$stunting_episode <- rep(NA, N)
  d$stunting_episode[d$stunt_inc==1] <- "stunted"
  d$stunting_episode[d$stunt_rec_inc==1] <- "Not stunted"
  
  
  #Have to mark first observations as stunted or not stunted if dropBornstunted=F
  if(dropBornstunted==F){
    d <- d %>% group_by(SUBJID) %>% 
      mutate(stunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & stunt==0, "Not stunted", stunting_episode),
             stunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & stunt==1, "stunted", stunting_episode),
             born_stunt_inc= 0,
             stunting_episode = na.locf(stunting_episode, fromLast=F)) %>% #Last observation carried forward 
      ungroup()
  }else{
    d <- d %>% group_by(SUBJID) %>% 
      mutate(stunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & stunt==0, "Not stunted", stunting_episode),
             stunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & stunt==1, "Born stunted", stunting_episode),
             stunt_inc = ifelse(stunting_episode=="Born stunted",0, stunt_inc),
             born_stunt_inc= ifelse(AGEDAYS==min(AGEDAYS) & stunting_episode=="Born stunted",1,0),
             stunting_episode = na.locf(stunting_episode, fromLast=F)) %>% #Last observation carried forward 
      ungroup()      
  }
  
  #Indicate risk of stunting or recovery 
  d$stunt_risk[(d$stunting_episode=="Not stunted" & d$past_stunt==0) | d$stunt_inc==1] <- 1 
  d$stunt_rec_risk[(d$stunting_episode!="Not stunted" & d$stunt_inc!=1) | d$stunt_rec_inc==1] <- 1 
  
  
  
  #Calculate duration of stunting episodes
  d <- d %>%  group_by(SUBJID) %>%
    mutate(episode_ID = cumsum(born_stunt_inc+stunt_inc+stunt_rec_inc) + 1) %>% #Create unique episode ID
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
  
  
  #Variable for duration of only stunting episodes
  d$stunting_duration <- NA
  d$stunting_duration[d$stunting_episode=="stunted"] <- d$duration[d$stunting_episode=="stunted"]
  
  
  
  #---------------------------------------------------------
  #Calculate severe stunting and severe stunting recovery incidence and risk
  #---------------------------------------------------------    
  
  #Mark severe stunting changes
  d$sevstunt_falter <- NA
  d$sevstunt_falter[d$sevstuntchange==1 & d$past_sevstunt==0] <- 1
  d$sevstunt_falter[d$sevstuntchange!=1 | d$past_sevstunt==1] <- 0
  
  d$born_sevstunt_inc<-0
  d$sevstunting_episode <- rep(NA, N)
  d$sevstunting_episode[d$sevstunt_falter==1] <- "Severe stunted"
  d$sevstunting_episode[d$stunt_rec_inc==1] <- "Not Severe stunted"
  
  #Have to mark first observations as stunted or not stunted if dropBornstunted=F
  if(dropBornstunted==T){
    d <- d %>% group_by(SUBJID) %>% 
      mutate(
        born_sevstunt_inc= ifelse(AGEDAYS==min(AGEDAYS) & sevstunt==1,1,0),
        sevstunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevstunt==0, "Not Severe stunted", sevstunting_episode),
        sevstunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevstunt==1, "Born Severe stunted", sevstunting_episode),
        sevstunting_episode = na.locf(sevstunting_episode, fromLast=F),
        sevstunting_episode_lag=lag(sevstunting_episode)) %>% #Last observation carried forward 
      ungroup()      
  }else{
    d <- d %>% group_by(SUBJID) %>% 
      mutate(
        sevstunting_episode = ifelse
        (AGEDAYS==min(AGEDAYS) & sevstunt==0, "Not Severe stunted", sevstunting_episode),
        sevstunting_episode = ifelse(AGEDAYS==min(AGEDAYS) & sevstunt==1, "Severe stunted", sevstunting_episode),
        sevstunting_episode = na.locf(sevstunting_episode, fromLast=F),
        sevstunting_episode_lag=lag(sevstunting_episode)) %>% #Last observation carried forward 
      ungroup()     
  }
  
  #Indicate incidence of severe stunting and recovery
  d$sevstunting_episode_lag[is.na(d$sevstunting_episode_lag)]<-"Not Severe stunted"
  d$sevstunt_inc <- d$sevstunt_rec_inc <- 0
  d$sevstunt_inc[d$sevstunting_episode=="Severe stunted" & (d$sevstunting_episode_lag=="Not Severe stunted")] <- 1
  d$sevstunt_rec_inc[d$sevstunting_episode=="Not Severe stunted" & (d$sevstunting_episode_lag!="Not Severe stunted")] <- 1
  
  
  #Create unique severe stunting episode IDs
  d <- d %>%  group_by(SUBJID) %>%
    mutate(sev_episode_ID = cumsum(sevstunt_inc+born_sevstunt_inc+sevstunt_rec_inc) + 1) %>% #Create unique episode ID
    ungroup() 
  
  #Fill in severe stunting  episodes
  # d <- d %>% group_by(SUBJID, sev_episode_ID) %>%
  #   mutate(sevstunting_episode = na.locf(sevstunting_episode, fromLast=F)) %>%
  #   ungroup()
  
  
  #Mark risk of severe stunting and severe stunting recovery 
  d$sevstunt_risk[(d$sevstunting_episode!="Severe stunted" & d$sevstunting_episode!="Born Severe stunted" & d$past_sevstunt==0) | d$sevstunt_inc==1] <- 1 
  d$sevstunt_rec_risk[((d$sevstunting_episode=="Severe stunted" | d$sevstunting_episode=="Born Severe stunted") & d$sevstunt_inc!=1) | d$sevstunt_rec_inc==1] <- 1 
  
  
  #Calculate duration of severe stunting episodes
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
  
  
  #Variable for duration of only severe stunting episodes
  d$sevstunting_duration <- NA
  if(dropBornstunted==F){
    d$sevstunting_duration[d$sevstunting_episode=="Severe stunted"] <- d$sevduration[d$sevstunting_episode=="Severe stunted"]
  }else{
    d$sevstunting_duration[d$sevstunting_episode=="Severe stunted" | d$sevstunting_episode=="Born Severe stunted"] <- d$sevduration[d$sevstunting_episode=="Severe stunted" | d$sevstunting_episode=="Born Severe stunted"]
  }
  
  
  #Calculate 30,60, 90 day recovery and faltering into severe stunting
  d$period_30d <- d$period_60d <- d$period_90d <- T
  d$stunt_rec30d <- d$stunt_rec60d <- d$stunt_rec90d <- NA
  d$stunt_rec30d[d$stunt_inc==1]  <- d$stunt_rec60d[d$stunt_inc==1] <- d$stunt_rec90d[d$stunt_inc==1] <- 0
  d$sevstunt_inc30d <- d$sevstunt_inc60d <- d$sevstunt_inc90d <- NA
  d$sevstunt_inc30d[d$stunt_inc==1] <- d$sevstunt_inc60d[d$stunt_inc==1] <- d$sevstunt_inc90d[d$stunt_inc==1] <- 0
  for(i in 1:90){
    d <- d %>% group_by(SUBJID) %>%
      mutate(
        rec_inc_lead_i = lead(stunt_rec_inc, i),
        sev_inc_lead_i = lead(sevstunt_inc, i),
        days_lead_i = abs(lead(midpoint_age, i) - midpoint_age)
      )
    
    d$period_30d[d$days_lead_i > 30] <- F
    d$period_60d[d$days_lead_i > 60] <- F
    d$period_90d[d$days_lead_i > 90] <- F
    
    d$stunt_rec30d[d$stunt_inc==1 & d$period_30d & d$rec_inc_lead_i==1] <- 1
    d$stunt_rec60d[d$stunt_inc==1 & d$period_60d & d$rec_inc_lead_i==1] <- 1
    d$stunt_rec90d[d$stunt_inc==1 & d$period_90d & d$rec_inc_lead_i==1] <- 1
    
    d$sevstunt_inc30d[d$stunt_inc==1 & d$period_30d & d$sev_inc_lead_i==1] <- 1
    d$sevstunt_inc60d[d$stunt_inc==1 & d$period_60d & d$sev_inc_lead_i==1] <- 1
    d$sevstunt_inc90d[d$stunt_inc==1 & d$period_90d & d$sev_inc_lead_i==1] <- 1
    
    d$stunt_rec30d[d$stunt_inc==1 & d$period_30d & d$rec_inc_lead_i==0 & d$stunt_rec30d!=1] <- 0
    d$stunt_rec60d[d$stunt_inc==1 & d$period_60d & d$rec_inc_lead_i==0 & d$stunt_rec60d!=1] <- 0
    d$stunt_rec90d[d$stunt_inc==1 & d$period_90d & d$rec_inc_lead_i==0 & d$stunt_rec90d!=1] <- 0
    
    d$sevstunt_inc30d[d$stunt_inc==1 & d$period_30d & d$sev_inc_lead_i==0 & d$sevstunt_inc30d!=1] <- 0
    d$sevstunt_inc60d[d$stunt_inc==1 & d$period_60d & d$sev_inc_lead_i==0 & d$sevstunt_inc60d!=1] <- 0
    d$sevstunt_inc90d[d$stunt_inc==1 & d$period_90d & d$sev_inc_lead_i==0 & d$sevstunt_inc90d!=1] <- 0
    
    #Stop for loop if all current leading observations are beyond 90 days
    if(min(d$days_lead_i, na.rm=T) > 90) break
  }
  
  #Drop intermediate variables
  d <- subset(d, select = -c(agelag, stuntlag, sevstuntlag, midpoint_age, stuntchange, sevstuntchange, past_stunt, past_sevstunt,
                             future_stunt, future_sevstunt,  sevstunt_falter, sevstunting_episode_lag, sev_incident_age, sev_maxage,
                             sevduration, rec_inc_lead_i, sev_inc_lead_i, days_lead_i, period_30d,period_60d,period_90d, next_midpoint
  )) %>%
    ungroup() %>% as.data.frame()
  if(dropBornstunted==T){
    d <- subset(d, select = -c(born_stunt_inc, born_sevstunt_inc)) %>%
      ungroup() %>% as.data.frame()      
  }
  
  #merge back in other columns
  d <- merge(d, othercolumns, by=c("SUBJID", "AGEDAYS"))
  
  #rename columns to match other functions
  d <- d %>% rename(stunt_rec = stunt_rec_inc,
                    sevstunt_rec = sevstunt_rec_inc) 
  
  return(d)
}  





#----------------------------------------------
#Function to calculate summary tables
#----------------------------------------------

stuntIncTable<-function(d, strat=T, agecats=c(6*30, 12*30, 18*30, 24*30), agecat_rownames=c("0-6 months","6-12 months", "12-18 months", "18-24 months")){
  
  
  if(strat==T){
    d$agecat <- as.factor(findInterval(d$AGEDAYS, agecats, rightmost.closed=F))
  }  
  
  summary<-stuntIncSummary(d, strat=strat)
  tab<-summary[[1]]
  means<-data.frame(strata=rep("Overall",nrow(summary[[2]])), summary[[2]])
  
  if(strat==T){
    
    strattab<-stratmeans<-NULL
    for(i in 1:(length(agecats))){
      temp<-stuntIncSummary(d[d$agecat==(i-1),], strat=F)
      strattab<-rbind(strattab, temp[[1]])
      stratmeans<-rbind(stratmeans, data.frame(strata=rep(agecat_rownames[i],nrow(temp[[2]])), temp[[2]]))
    }
    stratmeans<-stratmeans[!is.na(stratmeans$Mean),]
    
    #Replace truncated durations with stratified durations calculated from full data
    #numcats<-length(unique(stratmeans$strata))
    #if(numcats>nrow(summary[[3]])){numcats<-nrow(summary[[3]])}
    #strattab$total_duration[!is.na(strattab$total_duration)] <- summary[[3]]$total_duration[1:numcats]
    #strattab$average_duration[!is.na(strattab$total_duration)]  <- summary[[3]]$average_duration[1:numcats]
    #agelvls=levels(stratmeans$strat)
    #for(i in 1:(numcats)){
    #  stratmeans[stratmeans$strata== agelvls[i] & 
    #               stratmeans$statistic=="Average\nduration\nof\nstunting" ,3:8] <- summary[[4]][summary[[4]]$agecat==(i-1), -1]
    #}
    
    
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

stuntIncSummary<-function(d, strat=F){
  require(epitools)
  
  d <- d %>% arrange(SUBJID, AGEDAYS)
  
  #Average number of measurements per child
  child_nmeas <- d %>% group_by(SUBJID) %>%
    summarize(num_measurements=n()) %>% ungroup() %>% 
    summarize(num_children=n(),
              num_measurements=mean(num_measurements))
  
  
  #Overall sum and means
  Incidence_df <- d %>% 
    summarize(      
      prev_stunt= mean(HAZ < -2, na.rm=T) * 100,
      prev_sevstunt= mean(HAZ < -3, na.rm=T) * 100,
      anystunt= ifelse(sum(stunt) > 0 ,1,0) * 100,
      anysevstunt= ifelse(sum(sevstunt) > 0 ,1,0) * 100,
      persontime=sum(stunt_risk*period_length, na.rm=T),
      sev_persontime=sum(sevstunt_risk*period_length, na.rm=T),
      recovery_persontime=sum(stunt_rec_risk*period_length, na.rm=T),
      sevrecovery_persontime=sum(sevstunt_risk*period_length, na.rm=T),
      stunt_ep=sum(stunt_inc, na.rm=T),
      sevstunt_ep=sum(sevstunt_inc, na.rm=T),
      stunt_rec_ep=sum(stunt_rec, na.rm=T),
      sevstunt_rec_ep=sum(sevstunt_rec, na.rm=T),
      recoveries30d=sum(stunt_rec30d==1, na.rm=T),
      recoveries60d=sum(stunt_rec60d==1, na.rm=T),
      recoveries90d=sum(stunt_rec90d==1, na.rm=T),
      sevstunt30d=sum(sevstunt_inc30d==1, na.rm=T),
      sevstunt60d=sum(sevstunt_inc60d==1, na.rm=T),
      sevstunt90d=sum(sevstunt_inc90d==1, na.rm=T),
      no_recoveries30d=sum(stunt_rec30d==0, na.rm=T),
      no_recoveries60d=sum(stunt_rec60d==0, na.rm=T),
      no_recoveries90d=sum(stunt_rec90d==0, na.rm=T),
      no_sevstunt30d=sum(sevstunt_inc30d==0, na.rm=T),
      no_sevstunt60d=sum(sevstunt_inc60d==0, na.rm=T),
      no_sevstunt90d=sum(sevstunt_inc90d==0, na.rm=T)
    ) %>% 
    mutate(  stuntIR=stunt_ep/persontime * 1000,
             sevstuntIR=sevstunt_ep/sev_persontime * 1000,
             stuntrecIR=stunt_rec_ep/recovery_persontime * 1000,
             sevstuntrecIR=sevstunt_rec_ep/sevrecovery_persontime * 1000,
             perc_stuntrec_30d= sum(recoveries30d)/sum(no_recoveries30d)*100,
             perc_stuntrec_60d= sum(recoveries60d)/sum(no_recoveries60d)*100,
             perc_stuntrec_90d= sum(recoveries90d)/sum(no_recoveries90d)*100,
             perc_sevstuntinc_30d= sum(sevstunt30d)/sum(no_sevstunt30d)*100,
             perc_sevstuntinc_60d= sum(sevstunt60d)/sum(no_sevstunt60d)*100,
             perc_sevstuntinc_90d= sum(sevstunt90d)/sum(no_sevstunt90d)*100
    ) %>%
    as.data.frame()
  
  #Calculate average episode lengths
  #d <- d %>% group_by(SUBJID) %>%
  #  mutate(state_run = cumsum( stunt_inc+sevstunt_inc+stunt_rec) + 1,
  #         stunt_run = ifelse(stunt==1, state_run, 0)) %>%
  #  ungroup() %>%  group_by(SUBJID, state_run) %>% 
  #  mutate(state_dur = sum(period_length),
  #         stunt_dur = ifelse(stunt==1, state_dur, 0)) %>%
  #  as.data.frame()
  
  #data frame of episode durations
  #episode_duration <- d %>% 
  #  filter(stunt_inc==1) %>% #drop non-stunting periods
  #  group_by(SUBJID, state_run ) %>%
  #  slice(1) %>% ungroup() %>% 
  #  subset(., select=c(SUBJID,agecat,stunt_dur)) %>% 
  #  as.data.frame()
  
  episode_duration <- d %>% 
    filter(stunt_inc==1) %>% #drop non-stunting periods
    subset(., select=c(SUBJID,agecat,stunting_duration)) %>% 
    as.data.frame()
  
  # #Calculate mean, max, and total duration of stunting per child
  duration <- episode_duration %>%
    group_by(SUBJID) %>%
    summarize(
      total_duration= sum(stunting_duration, na.rm=T)
    )
  
  #Average episode length
  average_duration=mean(episode_duration$stunting_duration)
  #average total time stunted for each child
  total_duration=mean(duration$total_duration)
  
  #if(strat==T){
  #  duration_strat_average_duration <- episode_duration %>% group_by(agecat) %>% 
  #    summarize(average_duration=mean(stunt_dur))
  #  duration_strat_total_duration <- episode_duration %>% group_by(agecat, SUBJID) %>% 
  #    summarize(total_duration=sum(stunt_dur)) %>% 
  #    ungroup() %>% group_by(agecat) %>% 
  #    summarize(total_duration=mean(total_duration))
  #  duration_strat <- merge(duration_strat_average_duration, duration_strat_total_duration, by="agecat")
  #}
  
  
  #Create taable of summary statistics
  tab <- data.frame(child_nmeas, Incidence_df[1:4], total_duration, average_duration, Incidence_df[-c(1:4, 13:24)])
  
  #Calculate means
  means <- NULL
  try(
    means <- rbind(
      #longitudinal prevalence
      mean95CI(Y=(d$HAZ < -2), id=d$SUBJID, proportion=T, percent=F),
      mean95CI(Y=(d$HAZ < -3), id=d$SUBJID, proportion=T, percent=F),
      #duration
      mean95CI(Y=episode_duration$stunting_duration, id=episode_duration$SUBJID, proportion=F, percent=F),
      #incidence rates
      mean95CI(tab$stunt_ep, persontime=tab$persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$sevstunt_ep, persontime=tab$sev_persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$stunt_rec_ep, persontime=tab$recovery_persontime, count=T) * c(1,rep(1000,5)),
      mean95CI(tab$sevstunt_rec_ep, persontime=tab$sevrecovery_persontime, count=T) * c(1,rep(1000,5)),
      #percent recovery
      mean95CI(Y=d$stunt_rec30,  proportion=T, percent=F),
      mean95CI(Y=d$stunt_rec60,  proportion=T, percent=F),
      mean95CI(Y=d$stunt_rec90,  proportion=T, percent=F),
      mean95CI(Y=d$sevstunt_inc30d,  proportion=T, percent=F),
      mean95CI(Y=d$sevstunt_inc60d,  proportion=T, percent=F),
      mean95CI(Y=d$sevstunt_inc90d,  proportion=T, percent=F)))
  
  means <- data.frame(statistic =   c("Prevalence\nof\nstunting",
                                      "Prevalence\nof\nsevere\nstunting",
                                      "Average\nduration\nof\nstunting",
                                      "stunting\nincidence\nrate",
                                      "Severe\nstunting\nincidence\nrate",
                                      "stunting\nrecovery\nincidence\nrate",
                                      "Severe\nstunting\nrecovery\nincidence\nrate",
                                      "Percent\nstunting\nrecovered\nin 30 days",
                                      "Percent\nstunting\nrecovered\nin 60 days",
                                      "Percent\nstunting\nrecovered\nin 90 days",
                                      "Percent\nfalter to\nsevere\nstunting\nin 30 days",
                                      "Percent\nfalter to\nsevere\nstunting\nin 60 days",
                                      "Percent\nfalter to\nsevere\nstunting\nin 90 days"),means)
  
  #if(strat==T){
  #duration_strat_mean95CI <- NULL
  #for(i in levels(d$agecat)){
  # temp<-data.frame(agecat=i ,mean95CI(Y=episode_duration$stunt_dur[episode_duration$agecat==i], id=episode_duration$SUBJID[episode_duration$agecat==i]))
  # duration_strat_mean95CI <- rbind(duration_strat_mean95CI, temp)
  #}
  
  # return(list(tab, means, duration_strat, duration_strat_mean95CI))
  #}else{
  return(list(tab, means))
  #}
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
      prev_stunt,
      prev_sevstunt,
      anystunt,
      anysevstunt,
      #total_duration,
      average_duration
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Prevalence of stunting across all measurements` =  prev_stunt,
            `Prevalence of severe stunting across all measurements` = prev_sevstunt,
            `Proportion of children who were ever stunted`= anystunt,
            `Proportion of children who were ever severely stunted` = anysevstunt,
            #`Total duration of stunting episodes (days)` = total_duration,
            `Average duration of stunting episodes (days)` = average_duration)
  
  tab2 <- tab %>% 
    subset(., select=c(
      `Child age stratification`,
      num_children,
      num_measurements,
      stunt_ep,
      sevstunt_ep,
      persontime,
      sev_persontime,
      stuntIR,
      sevstuntIR
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Number of stunting episodes` = stunt_ep,
            `Number of severe stunting episodes` = sevstunt_ep,
            `No. of days at risk of stunting` = persontime,
            `No. of days at risk of severe stunting` = sev_persontime,
            `stunting incidence rate per 1000 days` = stuntIR,
            `Severe stunting incidence rate per 1000 days` = sevstuntIR
    )
  
  
  tab3 <- tab %>% 
    subset(., select=c(
      `Child age stratification`,
      num_children,
      num_measurements,
      stunt_rec_ep,
      recovery_persontime,
      stuntrecIR,
      perc_stuntrec_30d,       
      perc_stuntrec_60d,
      perc_stuntrec_90d
    )) %>%
    rename( `Number of children` = num_children,
            `Average number of measurements per child` = num_measurements,
            `Number of recoveries from stunting` = stunt_rec_ep,
            `No. of days at eligible for stunting recovery` = recovery_persontime,
            `stunting recovery incidence rate per 1000 days` = stuntrecIR,
            `Percent of stunting episodes recovered from in 30 days` = perc_stuntrec_30d,       
            `Percent of stunting episodes recovered from in 60 days` = perc_stuntrec_60d,
            `Percent of stunting episodes recovered from in 90 days` = perc_stuntrec_90d)
  return(list(tab1,tab2,tab3, tab))
}






#----------------------------------------------
#Plot functions
#----------------------------------------------

#GAM curve
HAZ_curve<-function(df){
  theme_set(theme_bw())
  
  # grab a color blind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(1,3,7)]
  
  p <- gmsn_age<-ggplot(df, aes(x = AGEDAYS)) +
    geom_smooth(aes(y=HAZ), color="#D55E00") +
    geom_jitter(aes(y=HAZ,x=AGEDAYS), height = 0.2, width=0.2,  alpha = 0.1, size=0.5)+
    labs(y = "HAZ",
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
  
  p <- ggplot(d, aes(x=AGEDAYS, y=HAZ)) + 
    geom_line() + guides(colour=FALSE) + xlab("Child Age (Days)") +
    ylab("HAZ") + aes(alpha=alpha, group=factor(SUBJID)) + guides(alpha=FALSE)
  
  return(p)
}


#Heatmap
heatmap<-function(d){
  theme_set(theme_bw())
  
  #generate HAZ categories
  d$stuntcat <- 0
  d$stuntcat[d$HAZ<(-2)] <- -1
  d$stuntcat[d$HAZ<(-3)] <- -2
  
  #Sum an ad-hoc total stunting score to rank children by  
  d <- d %>% group_by(SUBJID) %>% mutate(stuntscore=sum(stuntcat, na.rm=T)) %>% ungroup()
  table(d$stuntscore)
  
  #Add level for missingness
  d$stuntcat[is.na(d$HAZ)] <- 1
  d$stuntcat <- as.factor(d$stuntcat)
  d$SUBJID <- as.factor(d$SUBJID)
  table(d$stuntcat)
  
  #make ordered childnum by amount of stunting
  d <- d %>%
    arrange(stuntscore) %>%  
    mutate(SUBJID = factor(SUBJID, unique(SUBJID))) 
  d$childnum <- as.numeric(d$SUBJID)
  head(as.data.frame(d))
  
  #Create a child age in months variable
  d$agemonths<-floor(d$AGEDAYS/30)
  
  
  if(sum(d$stuntcat== -2) > 0){
    levels(d$stuntcat)<-c("Severely stunted", "stunted", "Not stunted")
    cbPalette <- c("#56B4E9", "#E69F00", "#c1c1c1")
  }else{
    levels(d$stuntcat)<-c( "stunted", "Not stunted")
    cbPalette <- c( "#E69F00", "#c1c1c1")   
  }
  
  p <- ggplot(d,aes(x=agemonths,y=childnum)) + 
    geom_tile(aes(fill=stuntcat)) + scale_fill_manual(values=cbPalette) +
    xlab("Child Age (Months)") + ylab("Child number")
  ggtitle("Heatmap of stunting and severe stunting episodes") +
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
# any stunting by 6 months and recovery in 60
# days among those stunted using TMLE, and "HR"/intensity
# ratio of incident stunting using pooled logistic regressions
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
                       born_not_stunt=F,
                       overall_dist=T,
                       run_anystunt=T,
                       run_recovery=T,
                       run_irr=T){
  
  cat("\n", unique(df$STUDYID),": ", unique(df$COUNTRY),"\n")
  
  Avar<-A
  
  #Rename risk factor to A
  colnames(df)[grep(paste0("^",A,"$"), colnames(df))]<-"A"
  
  #Drop any missing A
  df <- df %>% filter(!is.na(A))
  
  if(born_not_stunt==T){
    birth_df<-df %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
      slice(which.min(AGEDAYS)) %>% 
      mutate(bornstunt=as.numeric(HAZ < -2))
    birth_df$bornstunt[birth_df$AGEDAYS > 60] <- NA
    
    df <- left_join(df, birth_df, by= c("STUDYID", "SUBJID"))
    df <- df %>% filter(bornstunt==0)
  }
  
  #Tabulate n in agerange (to diagnose when no estimates are returned)
  n_in_agerange<- nrow(df[df$AGEDAYS>=agerange[1] & df$AGEDAYS<agerange[2],])
  
  inc <- df %>% group_by(SUBJID) %>% arrange(AGEDAYS) %>%
    mutate(stuntcount=cumsum(stunt_inc),
           anystunt=ifelse(stuntcount>0,1,0)) %>%
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
  
  if(run_anystunt==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),": ",dim(df),"\n")
    Anystunt6mo<-NULL
    try(
      Anystunt6mo<-tmle_risk(dat=as.data.frame(inc),
                            W=Wvars,
                            A="A",
                            n.cat=n.cat,
                            reflevel=reflevel,
                            Acuts=Acuts,
                            Alevels=Alevels,
                            outputdf=NULL,
                            Y="anystunt",
                            family="binomial",
                            SLlibrary=SLlibrary,
                            overall.dist=overall_dist,
                            sparseN=0,
                            adjusted=F))
    try(
      Anystunt6mo<-data.frame(study=rep(unique(df$STUDYID), nrow(Anystunt6mo)) ,country=rep(unique(df$COUNTRY), nrow(Anystunt6mo)) , variable=rep(Avar,  nrow(Anystunt6mo)),
                             Anystunt6mo, n_in_agerange=rep(n_in_agerange, nrow(Anystunt6mo))
                             , mean_age=rep(mean_age, nrow(Anystunt6mo))
                             , min_age=rep(min_age, nrow(Anystunt6mo))))
    
    if(is.null(Anystunt6mo)){
      Anystunt6mo<-data.frame(study=unique(df$STUDYID),country=unique(df$COUNTRY), variable=Avar,
                             variable.1=NA, level=NA,         ATE=NA,
                             ATE.var=NA,     ATE.CI1=NA,      ATE.CI2=NA,     ATE.Pval=NA,        RR=NA,    RR.CI1=NA,
                             RR.CI2=NA,      RR.Pval=NA,   logRR.psi=NA,   logRR.var=NA, compN=NA, refN=NA,    a=NA,    b=NA,
                             c=NA,    d=NA, meanLevel=NA, meanN=NA,     meanY=NA,   mean.sd=NA,     
                             mean.se=NA,  mean.CI1=NA, mean.CI2=NA, n_in_agerange=NA, mean_age=mean_age, min_age=min_age)
    }
  }else{
    Anystunt6mo<-NULL
  }
  
  
  if(run_recovery==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),"\n")
    drec<- df %>% filter(stunt_inc==1 & !is.na(A)  & AGEDAYS >= agerange[1]*30 & AGEDAYS < agerange[2]*30) %>%
      as.data.frame()
    
    if(Wvars[1]=="W1"){
      drec$W1<-rep(1, nrow(drec))
      drec$W2<-rep(1, nrow(drec))
    }
    
    #Need to add in SUBJID as clustering variable
    stuntrec60d<-NULL
    try(stuntrec60d<-tmle_risk(dat=as.data.frame(drec),
                              W=Wvars,
                              A="A",
                              n.cat=n.cat,
                              reflevel=reflevel,
                              Acuts=Acuts,
                              Alevels=Alevels,
                              outputdf=NULL,
                              Y="stunt_rec60d",
                              family="binomial",
                              SLlibrary="SL.glm",
                              overall.dist=overall_dist,
                              sparseN=0,
                              adjusted=F))
    try(
      stuntrec60d<-data.frame(study=rep(unique(df$STUDYID), nrow(stuntrec60d)) , country=rep(unique(df$COUNTRY), nrow(stuntrec60d)) , variable=rep(Avar,  nrow(stuntrec60d)),
                             stuntrec60d, n_in_agerange=rep(n_in_agerange, nrow(stuntrec60d))))
    
    if(is.null(stuntrec60d)){
      stuntrec60d<-data.frame(study=unique(df$STUDYID),country=unique(df$COUNTRY), variable=Avar,
                             variable.1=NA, level=NA,         ATE=NA,
                             ATE.var=NA,     ATE.CI1=NA,      ATE.CI2=NA,     ATE.Pval=NA,        RR=NA,    RR.CI1=NA,
                             RR.CI2=NA,      RR.Pval=NA,   logRR.psi=NA,   logRR.var=NA, compN=NA, refN=NA,    a=NA,    b=NA,
                             c=NA,    d=NA, meanLevel=NA, meanN=NA,     meanY=NA,   mean.sd=NA,     
                             mean.se=NA,  mean.CI1=NA, mean.CI2=NA, n_in_agerange=NA)
    }
    
  }else{
    stuntrec60d<-NULL
  }
  
  
  
  #pooled regression for stunting intensity
  if(run_irr==T){
    cat("\n",unique(as.character(df$STUDYID)),unique(as.character(df$COUNTRY)),"\n")
    stuntrisk <- df %>% filter(stunt_risk==1 & !is.na(SUBJID)  & !is.na(AGEDAYS)  & !is.na(stunt)  & AGEDAYS >= agerange[1]*30 & AGEDAYS < agerange[2]*30) %>% arrange(A)
    if(overall_dist==F){
      Acuts<-quantile(stuntrisk$A, probs = c((1:(n.cat-1))/n.cat), na.rm=T)
    }
    stuntrisk$Acat <- factor(findInterval(stuntrisk$A, Acuts, left.open=T))
    stuntrisk$Acat <- relevel(stuntrisk$Acat, ref=levels(stuntrisk$Acat)[reflevel])
    stuntrisk$ageZ<- scale(stuntrisk$AGEDAYS, center = TRUE, scale = TRUE)
    irr <- poissonRB(fmla="stunt~Acat + poly(ageZ,3,raw=TRUE)", dat=stuntrisk)
    logres<-irr$rfit
    
    res<-data.frame(study=rep(unique(df$STUDYID), nrow(logres)),country=rep(unique(df$COUNTRY), nrow(logres)) , variable=rep(Avar,  nrow(logres)), coef=rownames(logres), RR=exp(logres[,1]), ci.lb= exp(logres[,1]-1.96*logres[,2]), ci.ub=  exp(logres[,1]+1.96*logres[,2]),logRR=logres[,1] , log.se=logres[,2], pvalue=logres[,4])
  }else{
    res<-NULL
  }
  return(list(Anystunt6mo=Anystunt6mo, stuntrec60d=stuntrec60d, irr=res))
}












# stuntIncCalc_old<-function(d, washout=60){
#   require(tidyverse)
#   require(zoo)  
#   
#   
#   
#   #Filter out extreme or missing HAZ values
#   d <- d %>%  filter(!is.na(HAZ)) %>%
#     filter(HAZ > (-5) & HAZ < 5)
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
#   othercolumns <- d %>% subset(., select= -c(HAZ, VISITNUM)) 
#   d <- d %>% subset(., select= c(SUBJID, HAZ, AGEDAYS, VISITNUM)) 
#   
#   
#   #generate stunting
#   d$stunt= ifelse(d$HAZ < (-2),1,0)
#   d$sevstunt= ifelse(d$HAZ < (-3),1,0)
#   
#   d <- d %>%
#     arrange(SUBJID, AGEDAYS) %>%
#     group_by(SUBJID) %>%
#     mutate(stuntlag=lag(stunt),
#            sevstuntlag=lag(sevstunt),
#            agelag=lag(AGEDAYS),
#            agelead=lead(AGEDAYS),
#            period_lag= AGEDAYS - agelag, #Calculate number of days since last measurement
#            period_lead= agelead - AGEDAYS,
#            stuntchange= as.numeric((stunt - stuntlag)!=0 & !is.na(stuntlag)),
#            sevstuntchange=  as.numeric(((sevstunt - sevstuntlag)==1 & !is.na(sevstuntlag)))) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   #Fill in missing lag variables with 0 (for the first measurements)
#   d$stuntlag[is.na(d$stuntlag)] <- 0
#   d$sevstuntlag[is.na(d$sevstuntlag)] <- 0
#   d$agelag[is.na(d$agelag)] <- 0
#   d$period_lag[is.na(d$period_lag)] <- d$AGEDAYS[is.na(d$period_lag)] #Set days since birth
#   d$period_lead[is.na(d$period_lead)] <- 30  #Assume month if no further observations
#   d$period_length<-(d$period_lag + d$period_lead)/2 #Assumed period length
#   
#   #Expand longform dataset to have a row for every child day
#   vals <- expand.grid(AGEDAYS = 1:max(d$AGEDAYS),
#                       SUBJID = unique(d$SUBJID))
#   d <- left_join(vals, d, by=c("SUBJID", "AGEDAYS"))
#   dfill<-d %>% subset(.,select=c(AGEDAYS,SUBJID,stunt, sevstunt)) %>% 
#     group_by(SUBJID) %>%
#     na.locf(., fromLast=F) %>%
#     rename(stuntfill=stunt, sevstuntfill=sevstunt) %>% 
#     arrange(AGEDAYS) %>% ungroup()
#   d <- merge(d, dfill, by=c("SUBJID","AGEDAYS"))
#   d <- d %>% arrange(SUBJID,AGEDAYS)
#   
#   d <- d %>%
#     arrange(SUBJID,AGEDAYS) %>%
#     group_by(SUBJID) %>%
#     mutate(
#       daysstunt2m = roll_sum_fun(stuntfill, washout),
#       dayssevstunt2m = roll_sum_fun(sevstuntfill, washout),
#       pastdaysstunt2m = lag_sum_fun(stuntfill, washout),
#       pastdayssevstunt2m = lag_sum_fun(sevstuntfill, washout)) %>%
#     filter(!is.na(HAZ)) %>%
#     mutate( #Set past stunting to 0 for first observation
#       pastdaysstunt2m = ifelse(VISITNUM==1, 0, pastdaysstunt2m),
#       pastdayssevstunt2m = ifelse(VISITNUM==1, 0, pastdayssevstunt2m),
#       future_stunt = ifelse(daysstunt2m>0, 1,0),
#       future_sevstunt = ifelse(dayssevstunt2m>0, 1,0),
#       past_stunt = ifelse(pastdaysstunt2m>0, 1,0),
#       past_sevstunt = ifelse(pastdayssevstunt2m>0, 1,0)) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   #Calculate incidence
#   # Severe recovery incidence requires first full recovery from stunting
#   d <- d %>% group_by(SUBJID) %>%
#     mutate(
#       stunt_risk= as.numeric(stuntlag==0 & past_stunt==0),
#       sevstunt_risk= as.numeric(sevstuntlag==0 & past_sevstunt==0),
#       stunt_rec_risk= as.numeric(past_stunt==1),
#       sevstunt_rec_risk= as.numeric(past_sevstunt==1),
#       stunt_inc= as.numeric(stunt==1 & stuntlag==0 & past_stunt==0),
#       sevstunt_inc= as.numeric(sevstunt==1 & sevstuntlag==0 & past_sevstunt==0),
#       stunt_rec= as.numeric(stunt==0 & stuntlag==1 & future_stunt==0), 
#       times_stunted=cumsum(stunt_inc),
#       times_sevstunted=cumsum(sevstunt_inc),
#       times_rec=cumsum(stunt_rec),
#       times_rec_lag=lag(times_rec),
#       sevstunt_rec=ifelse(times_rec!=times_rec_lag & times_rec==times_sevstunted,1,0)
#     ) %>%
#     as.data.frame()
#   head(d,30)
#   
#   
#   
#   #Calculate 30,60, 90 day recovery
#   #And 30, 60, 90 day faltering into severe stunting
#   #Expand longform dataset to have a row for every child day
#   dsub <- d %>% subset(.,select=c(SUBJID, AGEDAYS, stunt_inc, stunt_rec, sevstunt_inc))
#   d_expand <- left_join(vals, dsub, by=c("SUBJID", "AGEDAYS"))
#   d_rec <- d_expand %>% group_by(SUBJID) %>% arrange(SUBJID,AGEDAYS) %>%
#     mutate(
#       stunt_rec30d= as.numeric(roll_sum_fun(stunt_rec, 30)>0 & (stunt_inc==1 & !is.na(stunt_inc))),
#       stunt_rec60d= as.numeric(roll_sum_fun(stunt_rec, 60)>0 & (stunt_inc==1 & !is.na(stunt_inc))),
#       stunt_rec90d= as.numeric(roll_sum_fun(stunt_rec, 90)>0 & (stunt_inc==1 & !is.na(stunt_inc))),     
#       
#       sevstunt_inc30d= as.numeric(roll_sum_fun(sevstunt_inc, 30)>0 & (stunt_inc==1 & !is.na(stunt_inc))),
#       sevstunt_inc60d= as.numeric(roll_sum_fun(sevstunt_inc, 60)>0 & (stunt_inc==1 & !is.na(stunt_inc))),
#       sevstunt_inc90d= as.numeric(roll_sum_fun(sevstunt_inc, 90)>0 & (stunt_inc==1 & !is.na(stunt_inc)))   
#     ) %>% 
#     filter(!is.na(stunt_inc)) %>%
#     subset(.,select=c(SUBJID, AGEDAYS, stunt_rec30d, stunt_rec60d, stunt_rec90d, sevstunt_inc30d, sevstunt_inc60d, sevstunt_inc90d)) %>%
#     ungroup() %>% as.data.frame()
#   
#   
#   d <- merge(d, d_rec, by=c("SUBJID","AGEDAYS"))
#   d <- d %>% arrange(SUBJID,AGEDAYS)
#   
#   #Drop intermediate variables
#   d <- d %>% subset(., select= -c(stuntlag, sevstuntlag, agelag, agelead, 
#                                   period_lag, period_lead, stuntchange, sevstuntchange, 
#                                   stuntfill, sevstuntfill, daysstunt2m, dayssevstunt2m, pastdaysstunt2m,
#                                   pastdayssevstunt2m, future_stunt, future_sevstunt, past_stunt, past_sevstunt))
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
# stuntIncSummary_old<-function(d, strat=F){
#   require(epitools)
#   
#   
#   #Calculated days at risk per person
#   Incidence_df <- d %>% group_by(SUBJID) %>%
#     summarize(persontime=sum(stunt_risk*period_length, na.rm=T),
#               sev_persontime=sum(sevstunt_risk*period_length, na.rm=T),
#               recovery_persontime=sum(stunt_rec_risk*period_length, na.rm=T),
#               sevrecovery_persontime=sum(sevstunt_risk*period_length, na.rm=T),
#               stunt_ep=sum(stunt_inc, na.rm=T),
#               sevstunt_ep=sum(sevstunt_inc, na.rm=T),
#               stunt_rec_ep=sum(stunt_rec, na.rm=T),
#               sevstunt_rec_ep=sum(sevstunt_rec, na.rm=T),
#               recoveries30d=sum(stunt_rec30d==1, na.rm=T),
#               recoveries60d=sum(stunt_rec60d==1, na.rm=T),
#               recoveries90d=sum(stunt_rec90d==1, na.rm=T),
#               sevstunt30d=sum(sevstunt_inc30d==1, na.rm=T),
#               sevstunt60d=sum(sevstunt_inc60d==1, na.rm=T),
#               sevstunt90d=sum(sevstunt_inc90d==1, na.rm=T),
#               no_recoveries30d=sum(stunt_rec30d==0, na.rm=T),
#               no_recoveries60d=sum(stunt_rec60d==0, na.rm=T),
#               no_recoveries90d=sum(stunt_rec90d==0, na.rm=T),
#               no_sevstunt30d=sum(sevstunt_inc30d==0, na.rm=T),
#               no_sevstunt60d=sum(sevstunt_inc60d==0, na.rm=T),
#               no_sevstunt90d=sum(sevstunt_inc90d==0, na.rm=T)
#     ) %>% 
#     mutate(  stuntIR=stunt_ep/persontime * 1000,
#              sevstuntIR=sevstunt_ep/sev_persontime * 1000,
#              stuntrecIR=stunt_rec_ep/recovery_persontime * 1000,
#              sevstuntrecIR=sevstunt_rec_ep/sevrecovery_persontime * 1000,
#              perc_stuntrec_30d= sum(recoveries30d)/sum(no_recoveries30d)*100,
#              perc_stuntrec_60d= sum(recoveries60d)/sum(no_recoveries60d)*100,
#              perc_stuntrec_90d= sum(recoveries90d)/sum(no_recoveries90d)*100,
#              perc_sevstuntinc_30d= sum(sevstunt30d)/sum(no_sevstunt30d)*100,
#              perc_sevstuntinc_60d= sum(sevstunt60d)/sum(no_sevstunt60d)*100,
#              perc_sevstuntinc_90d= sum(sevstunt90d)/sum(no_sevstunt90d)*100
#     ) %>%
#     as.data.frame()
#   #Replace perc_stuntrec with NA if no stunting occurred
#   Incidence_df$perc_stuntrec_30d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df$perc_stuntrec_60d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df$perc_stuntrec_90d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df$perc_sevstuntinc_30d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df$perc_sevstuntinc_60d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df$perc_sevstuntinc_90d[Incidence_df$stunt_ep==0]<-NA
#   Incidence_df
#   
#   
#   
#   #Calculate average episode lengths
#   d <- d %>% group_by(SUBJID) %>%
#     mutate(state_run = cumsum( stunt_inc+sevstunt_inc+stunt_rec) + 1,
#            stunt_run = ifelse(stunt==1, state_run, 0)) %>%
#     as.data.frame()
#   
#   d <- d %>% group_by(SUBJID, state_run) %>% 
#     mutate(state_dur = sum(period_length),
#            stunt_dur = ifelse(stunt==1, state_dur, 0)) %>%
#     as.data.frame()
#   
#   head(d,30)
#   
#   
#   #Calculate mean, max, and total duration of stunting per child
#   duration <- d %>% 
#     filter(stunt_run!=0) %>% #drop non-stunting periods
#     group_by(SUBJID, state_run ) %>% 
#     slice(1) %>% #Grab first row (so no duplication of period days)
#     ungroup() %>%
#     group_by(SUBJID) %>%
#     summarize(
#       total_duration= sum(stunt_dur, na.rm=T),
#       average_duration= mean(stunt_dur, na.rm=T))
#   head(duration)
#   
#   if(strat==T){
#     duration_strat <- d %>% 
#       filter(stunt_run!=0) %>% #drop non-stunting periods
#       group_by(SUBJID, state_run ) %>% 
#       slice(1) %>% #Grab first row (so no duplication of period days)
#       ungroup() %>%
#       group_by(SUBJID, agecat) %>%
#       summarize(
#         total_duration= sum(stunt_dur, na.rm=T),
#         average_duration= mean(stunt_dur, na.rm=T)) %>%
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
#       prev_stunt= mean(HAZ < -2, na.rm=T),
#       prev_sevstunt= mean(HAZ < -3, na.rm=T),
#       anystunt= ifelse(sum(stunt) > 0 ,1,0),
#       anysevstunt= ifelse(sum(sevstunt) > 0 ,1,0)
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
#   studymeans$stuntIR<-mean(sum(sumdf$stunt_ep, na.rm=T)/sum(sumdf$persontime), na.rm=T)*1000
#   studymeans$sevstuntIR<-mean(sum(sumdf$sevstunt_ep, na.rm=T)/sum(sumdf$sev_persontime), na.rm=T)*1000
#   studymeans$stuntrecIR<-mean(sum(sumdf$stunt_rec_ep, na.rm=T)/sum(sumdf$recovery_persontime), na.rm=T)*1000
#   studymeans$sevstuntrecIR<-mean(sum(sumdf$sevstunt_rec_ep, na.rm=T)/sum(sumdf$sevrecovery_persontime), na.rm=T)*1000
#   
#   
# 
#   tab<-cbind(nrow(sumdf) ,studymeans[,c(1:7,22:31)], studysum[,c(8:21)])
#   names(tab)[1]<-"num_children"
#   tab[3:6]<-tab[3:6]*100
# 
#   stuntIR.CI <- pois.exact(tab$stunt_ep, pt = tab$persontime, conf.level = 0.95)[3:5] * 1000
#   sevstuntIR.CI <- pois.exact(tab$sevstunt_ep, pt = tab$sev_persontime, conf.level = 0.95)[3:5] * 1000
#   stuntrecIR.CI <- pois.exact(tab$stunt_rec_ep, pt = tab$recovery_persontime, conf.level = 0.95)[3:5] * 1000
#   sevstuntrecIR.CI <- pois.exact(tab$sevstunt_rec_ep, pt = tab$sevrecovery_persontime, conf.level = 0.95)[3:5] * 1000
#   
#   mean_stunt_incIR <- data.frame(N=NA, Mean=stuntIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=stuntIR.CI[2] ,  Upper.95.CI=stuntIR.CI[3] )
#   mean_sevstunt_inc <- data.frame(N=NA, Mean=sevstuntIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=sevstuntIR.CI[2] ,  Upper.95.CI=sevstuntIR.CI[3] )
#   mean_stuntrecIR <- data.frame(N=NA, Mean=stuntrecIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=stuntrecIR.CI[2] ,  Upper.95.CI=stuntrecIR.CI[3] )
#   mean_sevstuntrecIR <- data.frame(N=NA, Mean=stuntrecIR.CI[1], SD=NA, Robust.SE=NA ,  Lower.95.CI=stuntrecIR.CI[2] ,  Upper.95.CI=stuntrecIR.CI[3] )
#   
#   colnames(mean_stunt_incIR) <- colnames(mean_sevstunt_inc) <- colnames(mean_stuntrecIR) <- colnames(mean_sevstuntrecIR) <- c("N","Mean","SD","Robust SE", "Lower 95%CI", "Upper 95%CI") 
#   
#   #Calculate means
#   means <- NULL
#   try(
#     means <- rbind(
#       mean95CI(Y=sumdf$prev_stunt, id=sumdf$SUBJID, proportion=T, percent=F),
#       mean95CI(Y=sumdf$prev_sevstunt, id=sumdf$SUBJID, proportion=T, percent=F),
#       mean95CI(Y=sumdf$average_duration, id=sumdf$SUBJID),
#       
#       mean_stunt_incIR,
#       mean_sevstunt_inc,
#       mean_stuntrecIR,
#       mean_sevstuntrecIR,
#       
#       mean95CI(Y=sumdf$perc_stuntrec_30d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_stuntrec_60d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_stuntrec_90d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevstuntinc_30d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevstuntinc_60d, id=sumdf$SUBJID, proportion=T, percent=T),
#       mean95CI(Y=sumdf$perc_sevstuntinc_90d, id=sumdf$SUBJID, proportion=T, percent=T)))
#   
#   means <- data.frame(statistic =   c("Prevalence\nof\nstunting",
#                                       "Prevalence\nof\nsevere\nstunting",
#                                       "Average\nduration\nof\nstunting",
#                                       "stunting\nincidence\nrate",
#                                       "Severe\nstunting\nincidence\nrate",
#                                       "stunting\nrecovery\nincidence\nrate",
#                                       "Severe\nstunting\nrecovery\nincidence\nrate",
#                                       "Percent\nstunting\nrecovered\nin 30 days",
#                                       "Percent\nstunting\nrecovered\nin 60 days",
#                                       "Percent\nstunting\nrecovered\nin 90 days",
#                                       "Percent\nfalter to\nsevere\nstunting\nin 30 days",
#                                       "Percent\nfalter to\nsevere\nstunting\nin 60 days",
#                                       "Percent\nfalter to\nsevere\nstunting\nin 90 days"),means)
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