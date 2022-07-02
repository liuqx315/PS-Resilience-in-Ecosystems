## Lable100	Area100	X100	Y100	Circularity100	EquivDiameter100	Orientation100	
# Lable110	Area110	X110	Y110	Circularity110	EquivDiameter110	Orientation110
remove(list=ls())
library(car)      # Anova, levene.test, outlier.test
library(nortest)  # ad.test, cvm.test
library(base)
require(lmerTest)     # For lmer()
require(car)          # For Anova (note the capital A)
library(Hmisc)
# library(lettercase)
library(sfsmisc)
library(gplots)
library(grDevices)
library(RColorBrewer)
library(stats)
library(MASS)
library(nloptr)
library(lme4)
library(broom)
library(emmeans)
library(tidyverse)
library(furrr)
library(ggplot2)
library(Rmisc)

stderr <- function(x) sqrt(var(x,na.rm=T)/length(which(!is.na(x))))
options(warn=-1)
# Define a function that opens a window on either Mac or Windows (Linux users use X11)
OpenWindow = function (Width,Height) {
  if (Sys.info()["sysname"]=="Darwin"){  # (Darwin stands for a Mac computer)
    quartz(width=Width, height=Height)
  } else {
    windows(width = Width, height = Height)}
}

SaveFigure=function(FileName){
  if (Sys.info()["sysname"]=="Darwin"){ # (Darwin stands for a Mac computer)
    quartz.save(paste(FileName,'.pdf',sep=''),type = c("pdf"),device = dev.cur())
  } else
    savePlot(filename = FileName,type = c("pdf"),device = dev.cur(),restoreConsole = TRUE) 
}

coul = brewer.pal(8, "Dark2") #Accent, Dark2,Paired,Pastel1,Pastel1,Set1,Set2,Set3,Spectral
mycolors <- brewer.pal(8, "Accent"); #terrain.colors(10)
#my_colors=colorRampPalette(my_colors)(100)
decol = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
          '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
          '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
          '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')
# Continuous colors
# scale_color_brewer(palette="Paired") 
# # Discrete colors
# scale_color_brewer(palette="Dark2") 
# # Gradient colors
# scale_color_brewer(palette="Accent") 

path=getwd()
setwd(path)


dataDDA = read.csv(file='DDA_100vs110.csv',header = T,sep = ",")
Tscale=10*20000/300 # time scale
# dataAll<-data.frame(dataDDA)
dat1<- dataDDA %>% mutate(GrowthArea=(Area110-Area100)/Tscale,
                         GrowthRadius=(EquivDiameter110-EquivDiameter100)/Tscale)
dat1<-subset(dat1,Area110>0)

# plot resutls
library(tidyverse)
library(modelr)
options(na.action = na.warn)
library(nycflights13)
library(lubridate)

pp<-ggplot(dat1,aes(y = GrowthArea, x = Area100)) + 
  geom_point(color= mycolors[1],shape=21, size = 2, alpha = 1)+
  geom_hline(yintercept= 0,linetype="dashed")+
  geom_smooth(method="loess", se=F,size=1) +
  # stat_smooth(method = "lm", formula = y ~ x*exp(-x), se=F, size = 1)+
  # scale_x_continuous(breaks = c(0,5,10,15,20),labels=c("0","5","10","15","20"))+
  # scale_color_manual(values = my_colors)+
  # scale_colour_manual(name = "Year",
  #                     labels = c("2013", "2014"),
  #                     values = c(coul[1], coul[2])) +   
  # scale_shape_manual(name = "Year",
  #                    labels = c("2013","2014"),
  #                    values = c(19,17))+
  scale_y_continuous(limits=c(-0.1,0.1)) +
  scale_alpha(guide=FALSE)+
  scale_size(guide=FALSE)+
  guides(color=FALSE)+
  guides(shape=FALSE)+
  theme_bw()+
  # guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  labs(x = "Patch size in area", y = "Patch growth rate", color = "")+
  theme(
    panel.border = element_rect(size = 1.0,colour = "black", linetype=1),
    # axis.line=element_line(size = 0.6, colour = "black", linetype=1),
    legend.position=c(0.86,0.85),
    legend.title = element_blank(), # element_blank()
    legend.text = element_text(size = 14,color = "black"),
    legend.background = element_rect(fill = NA,color = NA),
    axis.text.x = element_text(size=14, vjust = 0.5,hjust=0.5,color = 'black'),
    axis.text.y = element_text(size=14, hjust = 1,color = 'black'),
    axis.title.x = element_text(size=16,margin = margin(t = 5, b = 0)),
    axis.title.y = element_text(size=16,margin = margin(r = 2)),
    plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.12, "cm"))

print(pp)
ggsave("Fig4C.pdf", width = 5, height = 4)

## ============= plot SDF growth figure================

dataSDF = read.csv(file='SDF_50vs100.csv',header = T,sep = ",")
Tscale=50*3000/300 # time scale DeltT=50=100-50;
dat2<- dataSDF %>% mutate(GrowthArea=(Area110-Area100)/Tscale,
                         GrowthRadius=(EquivDiameter110-EquivDiameter100)/Tscale)
# dat<-subset(dat,Area110>0)

pp2<-ggplot(dat2,aes(y = GrowthArea, x = Area100)) + 
  geom_point(color= 'darkgrey',shape=21, size = 2, alpha = 1)+
  geom_hline(yintercept= 0,linetype="dashed")+
  geom_smooth(method="loess", se=F,size=1,  color='darkred') +
  # stat_smooth(method = "lm", formula = y ~ x*exp(-x), se=F, size = 1)+
  # scale_x_continuous(breaks = c(0,5,10,15,20),labels=c("0","5","10","15","20"))+
  # scale_color_manual(values = my_colors)+
  # scale_colour_manual(name = "Year",
  #                     labels = c("2013", "2014"),
  #                     values = c(coul[1], coul[2])) +   
  # scale_shape_manual(name = "Year",
  #                    labels = c("2013","2014"),
  #                    values = c(19,17))+
  scale_y_continuous(limits=c(-0.1,0.1)) +
  scale_alpha(guide=FALSE)+
  scale_size(guide=FALSE)+
  guides(color=FALSE)+
  guides(shape=FALSE)+
  theme_bw()+
  # guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  labs(x = "Patch size in area", y = "Patch growth rate", color = "")+
  theme(
    panel.border = element_rect(size = 1.0,colour = "black", linetype=1),
    # axis.line=element_line(size = 0.6, colour = "black", linetype=1),
    legend.position=c(0.86,0.85),
    legend.title = element_blank(), # element_blank()
    legend.text = element_text(size = 14,color = "black"),
    legend.background = element_rect(fill = NA,color = NA),
    axis.text.x = element_text(size=14, vjust = 0.5,hjust=0.5,color = 'black'),
    axis.text.y = element_text(size=14, hjust = 1,color = 'black'),
    axis.title.x = element_text(size=16,margin = margin(t = 5, b = 0)),
    axis.title.y = element_text(size=16,margin = margin(r = 2)),
    plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.12, "cm"))

print(pp2)
ggsave("Fig4D.pdf", width = 5, height = 4)
