remove(list=ls())
library(car)      # Anova, levene.test, outlier.test
library(nortest)  # ad.test, cvm.test
library(base)
require(lmerTest)     # For lmer()
require(car)          # For Anova (note the capital A)
# library(Hmisc)
# library(lettercase)
library(sfsmisc)
library(gplots)
library(grDevices)
library(RColorBrewer)
library(stats)
library(MASS)
library(nloptr)
library(lme4)
# library(broom)
library(emmeans)
# library(tidyverse)
library(furrr)
library(ggplot2)
library(Rmisc)
library(ggExtra)
library(scales)
require(mgcv)
# library(modelr)
# options(na.action = na.warn)
library(nycflights13)
library(lubridate)
library(splines)
library(MASS)
library(data.table)
library(tikzDevice)

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

revlog_trans <- function(base = exp(1)){
  ## Define the desired transformation.
  trans <- function(x){
    -log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(-x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}


dataAll1 = read.csv(file='DDA_LS.csv',header = F,sep = "\t")
dataAll2 = read.csv(file='SDF_LS.csv',header = F,sep = "\t")

datMean=data.frame(cbind(t(apply(dataAll1[-1,], 2, mean)),t(apply(dataAll2[-1,], 2, mean))))
datSD=data.frame(cbind(t(apply(dataAll1[-1,], 2, sd)),t(apply(dataAll2[-1,], 2, sd))))
Time=data.frame(cbind(dataAll1[1,],dataAll2[1,]))

Model<-rep(c("DDA model","SDFs model"),c(dim(dataAll1)[2],dim(dataAll2)[2]))

data=data.frame(transpose(Time),transpose(datMean),transpose(datSD),Model)
colnames(data)=c('Time','datMean','datSD','Model')

# plot resutls
# fun.0 <- function(x) exp(1.56-0.29*x)
# fun.1 <- function(x) 1.7*(x-2)*exp(-0.3*(x-2))-0.5
# fun.2 <- function(x) 2.2*(x-2)*exp(-0.2*x)*(1-x/12)
fun.3 <- function(x) 5*x^{0.28}

pp<-ggplot(data,aes(y=datMean, x=Time,color= Model,shape= Model)) +
  geom_line(aes(y=datMean, x=Time,color= Model, shape= Model), size = 1, alpha = 1)+
  geom_hline(yintercept= 17,linetype="dashed")+
  geom_smooth(method="loess", se=T, size=1.0,fill = "grey50",span = 0.05,
              alpha=0.3, linetype = 1) + # "lm", "glm", "gam", "loess"
  # geom_point(aes(color= Model, shape= Model), size = 2, alpha = 1)+
  # annotate("segment", x = 5, xend = 30, y = 3e-7, yend = 4e-4, colour = "purple", 
  #          size=1, alpha=1.0,linetype="dashed")+
  geom_ribbon(data=data,aes(ymin=datMean-1.0*datSD,
                            ymax=datMean+1.0*datSD),size=0.0,alpha=0.2)+
  stat_function(fun = fun.3,size=1.0,col='steelblue',linetype = 2,aes(color = "Fitted") ) +
  # scale_x_continuous(breaks = c(0,5,10,15,20),labels=c("0","5","10","15","20"))+
  scale_x_continuous(trans = "log10",breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)),limits=c(8,1e5)) +
  scale_y_continuous(trans = 'log10',breaks = c(10,20,40,100),limits=c(10,120)) +
  # scale_x_continuous(breaks = c(1,10,100,1000), trans = revlog_trans(base = 10))+
  scale_fill_brewer(palette = "Accent")+
  # scale_shape_manual(name = "tau",
  #                    labels = c("3.0","4.5","5.0"),
  #                    values = c(1,2,0))+
  scale_alpha(guide=FALSE)+ scale_size(guide=FALSE)+ guides(shape=FALSE)+ guides(fill=FALSE)+
  # guides(color=FALSE)+
  # guides(color=guide_legend(nrow = 3, byrow = T,override.aes = list(size=1.0)))+
  # guides(shape = guide_legend(nrow = 3, byrow = T,
  #                             override.aes = list(size=2.5)))+
  annotation_logticks(sides="lb",short = unit(0.10, "cm"),  mid = unit(0.15, "cm"),
                      long = unit(0.20, "cm"),size = 0.3)+
  theme_bw()+
  labs(x="time", y='wavelength', color = "")+
  theme(
    panel.border = element_rect(size = 1.0,colour = "black", linetype=1),
    # axis.line=element_line(size = 0.6, colour = "black", linetype=1),
    legend.position=c(0.25,0.9),
    legend.title = element_blank(), # element_blank()
    legend.text = element_text(size = 12,color = "black"),
    legend.background = element_rect(fill = NA,color = NA),
    legend.spacing.y = unit(0., 'cm'),
    axis.text.x = element_text(size=14, vjust = 0.5,hjust=0.5,color = 'black'),
    axis.text.y = element_text(size=14, hjust = 1,color = 'black'),
    axis.title.x = element_text(size=16,margin = margin(t = 5, b = 0)),
    axis.title.y = element_text(size=16,margin = margin(r = 2)),
    axis.ticks.length = unit(.0, "cm"),
    plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.12, "cm"))

print(pp)

ggsave("FigS12DLS.pdf", width = 5, height = 4)

