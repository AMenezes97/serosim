## Load necessary packages
library(ggplot2)
library(tidyverse)
library(cowplot)

## Import run time data set 
run_time_path <- system.file("run_time", "serosim_run_times_new.csv", package = "serosim")
df <- read.csv(file = run_time_path, header = TRUE)

## Separate run time by case study 
df_readme <-  df %>% filter(df$example=="README")
df_cs1<-df %>% filter(df$example=="CS1")
df_cs2 <-df %>% filter(df$example=="CS2")

theme_set(theme_bw()) 

## Create each individual plot 
p_r <- ggplot(df_readme, aes(x=as.factor(individuals), y=(mean/60), colour=as.character(times)))
p_r <- p_r + geom_point(size=2)  +
  geom_errorbar(aes(ymin=range_min/60,ymax=range_max/60)) +
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
  guides(color = guide_legend(title = "Number of time steps"))

p_1 <- ggplot(df_cs1, aes(x=as.factor(individuals), y=(mean/60), colour=as.character(times)))
p_1 <- p_1+ geom_point(size=2)  +
  geom_errorbar(aes(ymin=range_min/60,ymax=range_max/60)) +
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
  guides(color = guide_legend(title = "Number of time steps"))

p_2 <- ggplot(df_cs2, aes(x=as.factor(individuals), y=(mean/60), colour=as.character(times)))
p_2 <- p_2+ geom_point(size=2)  +
  geom_errorbar(aes(ymin=range_min/60,ymax=range_max/60)) +
  labs(y="Run time in minutes", 
       x="Number of individuals", 
       title="Case study 2 ",
       subtitle = "3 exposure events, 2 biomarkers",
       key ="test") +
  theme(plot.title = element_text(hjust = 0.5, size=17)) +
  theme(axis.text.x = element_text(vjust=0.6, size= 12)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 12)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(plot.subtitle = element_text(hjust=0.5, size= 12)) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(title = "Number of time steps"))

## Combine all 3 plots
## Export as pdf with width of 15 inches 
plot_grid(p_r,p_1,p_2, nrow=1, ncol=3, align = "hv", scale=c(.98,.98, .98))

