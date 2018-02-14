
#-----------------------------------
# wasting-study-inventory-heatmaps.R
#
# andrew mertens (amertens@berkeley.edu)
#
# create a heatmap of data availability
# for weight and height by Indian
# studies on GHAP
#
#-----------------------------------

#-----------------------------------
# input files:
#    heatmap_df.Rdata
#
# output files:
#    
#-----------------------------------

#-----------------------------------
# preamble
#-----------------------------------
rm(list=ls())
library('tidyverse')
library('stringr')
library('scales')
library('RColorBrewer')
library('gridExtra')


# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"


#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
  "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

#define a color for fonts
textcol <- "grey20"

#Plot theme
theme_set(theme_bw())


#-----------------------------------
# load the monthly summary data
#-----------------------------------
load('C:/Users/andre/Documents/HBGDki/Results/Rally 6A/heatmap_df.Rdata')
load('C:/Users/andre/Documents/HBGDki/Results/Rally 6A/metadata6A.Rdata')

#Drop studies without anthro or cross-sectional studies from the metadata
meta$study_id[!is.na(meta$notes)]

meta <- meta$study_id[is.na(meta$notes)]


#Drop non-India cohorts from multi-site studies:
df <-df[df$COUNTRY=="INDIA",]

#drop obs. of children over 4 years old
df <- df[as.numeric(as.character(df$agemonth))<48,]

#Tablulate number of months and overall mean anthro for each study
tb<-df %>% group_by(STUDYID) %>% summarize(Nmonths=n(), aveWHZ=mean(meanWHZ, na.rm=T))

df <- left_join(df, tb, by="STUDYID")



#Drop cross-sectional studies
df$STUDYID[df$Nmonths==1]

df <- df[df$Nmonths>1,]

#Drop studies without anthropometry
df$STUDYID[is.na(df$aveWHZ)]

df <- df[!(is.na(df$aveWHZ)),]




#Order studies by age range of children
df <- df %>% arrange(-Nmonths)

#Drop missing study
df <- df %>% filter(!is.na(STUDYID))

#remove grant identifiers
df$STUDYID<- gsub("^k.*?-" , "", df$STUDYID)
df$STUDYID <- factor(df$STUDYID, levels=unique(df$STUDYID))

#Make N categories
summary(df$N)
df$Ncat <- cut(df$N,breaks=c(0,10,50,150,250,500,1000,45574),labels=c("<10","10-50","50-150","150-250","250-500","500-1000",">1000"))
df$Ncat <- factor(df$Ncat)

#Set mean anthro to NA when less than 10 observations
df$meanWHZ[df$Ncat=="<10"] <- NA
df$meanHAZ[df$Ncat=="<10"] <- NA
df$wast[df$Ncat=="<10"] <- NA
df$stunt[df$Ncat=="<10"] <- NA


# heat map plot scheme
hm <- ggplot(df,aes(x=agemonth, y=STUDYID)) +
  # facet over measurement frequency
  geom_tile(colour="white",size=0.25)+
  #remove extra space
  #scale_y_discrete(expand=c(0,0))+
  # scale_x_continuous(expand=c(0,0),
  #                    breaks=1:7,labels=1:7)+
  #one unit on x-axis is equal to one unit on y-axis.
  #equal aspect ratio x and y axis
  coord_equal()+
  #set base size for all font elements
  theme_bw() +
  #theme options
  theme(
    # legend options
    legend.title=element_text(color=textcol,size=8),
    #reduce/remove legend margin
    legend.margin = margin(grid::unit(0.1,"cm")),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    legend.key.height=grid::unit(0.2,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(1,"cm"),
    #move legend to the bottom
    legend.position = "bottom",
    #set x axis text size and colour
    axis.text.x=element_text(size=8,colour=textcol,angle=0,vjust=0.5),
    #set y axis text colour and adjust vertical justification
    axis.text.y=element_text(size=8,vjust = 0.2,colour=textcol),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    # axis.ticks.x=element_blank(),
    #change title font, size, colour and justification
    plot.title=element_text(colour=textcol,hjust=0,size=12,face="bold"),
    #format facet labels
    strip.text.x = element_text(size=10),
    strip.text.y = element_text(angle=270,size=10),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank()
    #remove plot margins
    # plot.margin=margin(grid::unit(1,"cm"))
  )




N_hm <- hm +
  aes(fill=Ncat) +
  labs(x="Child age in months",y="",title="Number observations by study and child age") +
   scale_fill_brewer(palette = "Greens",na.value="grey90",
                    guide=guide_legend(title="N obs.",title.vjust = 1,
                                       label.position="bottom",label.hjust=0.5,nrow=1))
 
N_hm




whzhm <- hm +
  aes(fill=meanWHZ) +
  labs(x="Child age in months",y="",title="Mean WHZ ") +
  scale_fill_gradient(low = tableau10[2], high = tableau10[3],na.value="grey90",
                    guide=guide_legend(title="WHZ",title.vjust = 1,
                                       label.position="bottom",label.hjust=0.5,nrow=1))
 
whzhm


hazhm <- hm +
  aes(fill=meanHAZ) +
  labs(x="Child age in months",y="",title="Mean HAZ ") +
  scale_fill_gradient(low = tableau10[2], high = tableau10[10],na.value="grey90",
                    guide=guide_legend(title="HAZ",title.vjust = 1,
                                       label.position="bottom",label.hjust=0.5,nrow=1))
 
hazhm



#CATEGORIZE
wasthm <- hm +
  aes(fill=wast) +
  labs(x="Child age in months",y="",title="Wasting Prevalence") +
  scale_fill_gradient(low = tableau10[3], high = tableau10[2],na.value="grey90",
                    guide=guide_legend(title="Prevalence",title.vjust = 1,
                                       label.position="bottom",label.hjust=0.5,nrow=1))
 
wasthm


stunthm <- hm +
  aes(fill=stunt) +
  labs(x="Child age in months",y="",title="Stunting Prevalence ") +
  scale_fill_gradient(low = tableau10[10], high = tableau10[2],na.value="grey90",
                    guide=guide_legend(title="Prevalence",title.vjust = 1,
                                       label.position="bottom",label.hjust=0.5,nrow=1))
 
stunthm


#Creat 4 panel figure with Multiplot















