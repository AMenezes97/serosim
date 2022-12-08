## Run case study 1 (cs1_vignette.Rmd)

## Pull necessary simulation data 
ifxn_exposed<- res$exposure_histories_long %>%
  filter(value==1 & x==1) %>% 
  mutate(ifxn_time=t) %>% 
  select(i,ifxn_time)

vacc_exposed <- res$exposure_histories_long %>%
  filter(value==1 & x==2)  %>% 
  mutate(vacc_time=t) %>% 
  select(i,vacc_time)

obs60<-res$observed_antibody_states %>% 
  filter(t==60) %>% 
  mutate(obs_60=observed) %>% 
  select(i,obs_60)

obs120<-res$observed_antibody_states %>% 
  filter(t==120) %>% 
  mutate(obs_120=observed) %>% 
  select(i,obs_120)

## Combine all three datasets 
obs<- left_join(obs60,obs120,by="i")
df<- left_join(obs,ifxn_exposed,by="i")
df<- left_join(df,vacc_exposed,by="i")

## Create a seroconversion and seropositivity variable 
full_df <- df %>% 
  mutate(titer_diff=obs_120-obs_60) %>% 
  mutate(seroconv=ifelse(ifxn_time>=61 | vacc_time >=61,"yes","no")) %>% 
  mutate(seropos=ifelse(!is.na(ifxn_time) | !is.na(vacc_time),"yes","no")) 
save<-full_df


## Create data frame with thresholds to try 
thresholds <- seq(from=1, to=1000,by=.1)
sens_spec<- data.frame(
  threshold=thresholds,
  sensit=NA,
  specif=NA
)

for(thds in 1:nrow(sens_spec)){
  ## Pull the titer threshold for this analysis
  titer<- sens_spec$threshold[thds]
  
  ## Add a column to record whether each individual will be classified as seropositive or seronegative given the titer threshold
  data<-full_df %>% 
    mutate(test=ifelse(obs_120>=titer,"yes","no"))
  
  ## Add a new column to record whether each individual is a true/false positive or a true/false negative
  new_data<- data %>% 
    mutate(outcome=ifelse(test=="yes" & seropos=="yes","true positive",
                          ifelse(test=="yes" & seropos=="no", "false positive",
                                 ifelse(test=="no" & seropos=="yes","false negative",
                                        ifelse(test=="no" & seropos=="no","true negative",NA)))))
  ## Sum up the results 
  truepos<- sum(new_data$outcome=="true positive")
  falsepos<- sum(new_data$outcome=="false positive")
  trueneg<- sum(new_data$outcome=="true negative")
  falseneg<- sum(new_data$outcome=="false negative")
 
  ## Calculate sensitivity and specificity
  sens<- truepos/(truepos+falseneg)
  spec <- trueneg/(trueneg+falsepos)
  
  ## Save values in table
  sens_spec$sensit[thds]<-sens*100
  sens_spec$specif[thds]<-spec*100
}


## Create plot for paper
library(ggplot2)
library(lubridate)
theme_set(theme_bw())

p1<-ggplot(sens_spec, aes(x=threshold)) + 
  geom_line(aes(y=sensit, col="Sensitivity")) + 
  geom_line(aes(y=specif, col="Specificity")) + 
  labs(title="Sensitivity and specificity for varying thresholds of seropositivity", 
       subtitle="Case study 1: Measles", x="Threshold of Seropositivity (mIU/mL)", y="Percentage") +  # title and caption
  scale_color_manual(name="", 
                     values = c("Sensitivity"="#00ba38", "Specificity"="#f8766d")) +  # line color
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
  ggplot2::theme(plot.subtitle = element_text(hjust = 0.5, size=12)) +
  ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(legend.position = c(0.84, 0.5)) +
  theme(legend.text = element_text(colour="black", size=15)) + theme(legend.background=element_blank())
  



## Examine the sensitivity and specificity for setting the threshold of seroconversion 

## Create data frame with thresholds to try 
serocon_thresholds <- seq(from=1, to=300,by=.1)
sens_spec_serocon<- data.frame(
  serocon_threshold=serocon_thresholds,
  sensit=NA,
  specif=NA
)

for(thds in 1:nrow(sens_spec_serocon)){
  ## Pull the titer serocon.thresholds for this analysis
  titer_diff_thds<- sens_spec_serocon$serocon_threshold[thds]
  
  ## Add a column to record whether each individual will be classified as seropositive or seronegative given the titer serocon.thresholds
  data<-full_df %>% 
    mutate(test=ifelse(titer_diff>=titer_diff_thds,"yes","no"))
  
  ## Add a new column to record whether each individual is a true/false positive or a true/false negative
  new_data<- data %>% 
    mutate(outcome=ifelse(test=="yes" & seropos=="yes","true positive",
                          ifelse(test=="yes" & seropos=="no", "false positive",
                                 ifelse(test=="no" & seropos=="yes","false negative",
                                        ifelse(test=="no" & seropos=="no","true negative",NA)))))
  ## Sum up the results 
  truepos<- sum(new_data$outcome=="true positive")
  falsepos<- sum(new_data$outcome=="false positive")
  trueneg<- sum(new_data$outcome=="true negative")
  falseneg<- sum(new_data$outcome=="false negative")
  
  ## Calculate sensitivity and specificity
  sens<- truepos/(truepos+falseneg)
  spec <- trueneg/(trueneg+falsepos)
  
  ## Save values in table
  sens_spec_serocon$sensit[thds]<-sens*100
  sens_spec_serocon$specif[thds]<-spec*100
}


## Create plot for paper
library(ggplot2)
library(lubridate)
theme_set(theme_bw())

p2<-ggplot(sens_spec_serocon, aes(x=serocon_threshold)) + 
  geom_line(aes(y=sensit, col="Sensitivity")) + 
  geom_line(aes(y=specif, col="Specificity")) + 
  labs(title="Sensitivity and specificity for varying thresholds of seroconversion", 
       subtitle="Case study 1: Measles", x="Seroconversion Threshold (mIU/mL)", y="Percentage") +  # title and caption
  scale_color_manual(name="", 
                     values = c("Sensitivity"="#00ba38", "Specificity"="#f8766d")) +  # line color
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5, size=15)) +
  ggplot2::theme(plot.subtitle = element_text(hjust = 0.5, size=12)) +
  ggplot2::theme(axis.text.x = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.text.y = element_text(vjust=0.6, size= 10)) +
  ggplot2::theme(axis.title.y = element_text(vjust=0.6, size= 13)) +
  ggplot2::theme(axis.title.x = element_text(vjust=0.6, size= 13)) +
  theme(legend.position = c(0.8, 0.6)) +
  theme(legend.text = element_text(colour="black", size=15)) + theme(legend.background=element_blank())

plot_grid(p1, p2, nrow=1, ncol=2, align = "hv", scale=c(.98,.98))
