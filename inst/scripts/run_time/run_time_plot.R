## Load necessary packages
library(ggplot2)
library(tidyverse)
library(cowplot)

## Import run time data set 
# run_time_path <- system.file("run_time", "serosim_run_time_updated.csv", package = "serosim")
# df <- read.csv(file = run_time_path, header = TRUE)

library(readr)
serosim_run_time_updated <- read_csv("inst/scripts/run_time/serosim_run_time_updated.csv")
df<-serosim_run_time_updated


## If run times are over 120 seconds or are NA then set to 120 seconds
## Times that are NA are because the initial runs were over 120 seconds and so the 100 simulations were not done
df$min<-ifelse(df$min>120 | is.na(df$min),120,df$min)
df$max<-ifelse(df$max>120 | is.na(df$max),120,df$max)
df$median<-ifelse(df$median>120| is.na(df$median),120,df$median)
df$mean<-ifelse(df$mean>120| is.na(df$mean),120,df$mean)
df$first_q<-ifelse(df$first_q>120| is.na(df$first_q),120,df$first_q)
df$third_q<-ifelse(df$third_q>120| is.na(df$third_q),120,df$third_q)




## Separate run time by case study 
df_readme <-  df %>% filter(df$example=="quick_start")
df_cs1<-df %>% filter(df$example=="CS1")
df_cs2 <-df %>% filter(df$example=="CS2")
df_cs3 <-df %>% filter(df$example=="CS3")

theme_set(theme_bw()) 

## Create each individual plot 
p_r<- ggplot(df_readme, aes(x=as.factor(individuals), y=(mean/60),  colour=as.factor(times)))  
p_r<- p_r + geom_boxplot(aes(x=as.factor(individuals), ymin=(min/60), ymax=(max/60), lower= (first_q/60), middle=(median/60), upper=(third_q/60), group=interaction(as.factor(individuals), as.factor(times))), stat = "identity", width=0.5) +
  labs(y="Run time in minutes", 
       x="Number of individuals", 
       title="README example",
       subtitle = "2 exposure events, 1 biomarker",
       key ="test") +
  theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.text.x = element_text(vjust=0.6, size= 12)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 12)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(plot.subtitle = element_text(hjust=0.5, size= 12)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(title = "Number of time steps")) + scale_color_viridis_d(option="D", end=0.65)+
  scale_y_continuous(labels = c(0,0.5,1,1.5,">2"))




p_1<- ggplot(df_cs1, aes(x=as.factor(individuals), y=(mean/60),  colour=as.factor(times)))  
p_1<- p_1 +geom_boxplot(aes(x=as.factor(individuals), ymin=(min/60), ymax=(max/60), lower= (first_q/60), middle=(median/60), upper=(third_q/60), group=interaction(as.factor(individuals), as.factor(times))), stat = "identity", width=0.5) +
  labs(y="Run time in minutes", 
       x="Number of individuals", 
       title="Case study 1",
       subtitle = "2 exposure events, 1 biomarker",
       key ="test") +
  theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.text.x = element_text(vjust=0.6, size= 12)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 12)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(plot.subtitle = element_text(hjust=0.5, size= 12)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(title = "Number of time steps")) + scale_color_viridis_d(option="D", end=0.65)+
  scale_y_continuous(labels = c(0,0.5,1,1.5,">2"))


p_2<- ggplot(df_cs2, aes(x=as.factor(individuals), y=(mean/60),  colour=as.factor(times)))  
p_2<- p_2 +geom_boxplot(aes(x=as.factor(individuals), ymin=(min/60), ymax=(max/60), lower= (first_q/60), middle=(median/60), upper=(third_q/60), group=interaction(as.factor(individuals), as.factor(times))), stat = "identity", width=0.5) +
  labs(y="Run time in minutes", 
       x="Number of individuals", 
       title="Case study 2",
       subtitle = "3 exposure events, 2 biomarkers",
       key ="test") +
  theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.text.x = element_text(vjust=0.6, size= 12)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 12)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(plot.subtitle = element_text(hjust=0.5, size= 12)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(title = "Number of time steps")) + scale_color_viridis_d(option="D", end=0.65)+
  scale_y_continuous(labels = c(0,0.5,1,1.5,">2"))


## Create additional plot for case study 3
p_3<- ggplot(df_cs3, aes(x=as.factor(individuals), y=(mean/60),  colour=as.factor(times)))
p_3<- p_3 +geom_boxplot(aes(x=as.factor(individuals), ymin=(min/60), ymax=(max/60), lower= (first_q/60), middle=(median/60), upper=(third_q/60), group=interaction(as.factor(individuals), as.factor(times))), stat = "identity", width=0.5) +
  labs(y="Run time in minutes", 
       x="Number of individuals", 
       title="Case study 3",
       subtitle = "10 exposure events, 10 biomarkers",
       key ="test") +
  theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.text.x = element_text(vjust=0.6, size= 12)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 12)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(plot.subtitle = element_text(hjust=0.5, size= 12)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(title = "Number of time steps")) + scale_color_viridis_d(option="D", end=0.65)+
  scale_y_continuous(labels = c(0,0.5,1,1.5,">2"))


## Combine all 4 plots
plot_grid(p_r,p_1,p_2,p_3, nrow=2, ncol=2, align = "hv", scale=c(.98,.98, .98, .98))
## Export 10 x 10 
