setwd("~/Documents/GitHub/serosim")
library(tidyverse)
devtools::load_all()
devtools::document()

## Let's simulate 1 pathogen (pertussis) with the Teunis-style antibody model
## We'll generate the FOI from an SIR model


## 1. Simulation settings
times <- seq(1,365*2,by=1) ## 2 years in daily resolution
simulation_settings <- list("t_start"=min(times),"t_end"=max(times))

## 2. Population demography
N <- 1000
birth_times <- simulate_birth_times(N, times, age_min=0)
removal_times <- simulate_removal_times(N, times=times,birth_times=birth_times, removal_min=365,removal_max=max(times),prob_removal=1)

aux <- list("Sex"=list("name"="sex","options"=c("male", "female"), "proportion"=c(0.5,0.5)),
            "Group"=list("name"="group","options"=c("1", "2", "3", "4"), "proportion"=c(0.25,0.25,0.25,0.25)),
            "SES"=list("name"="ses","options"=c("low","high"),"proportion"=c(0.8,0.2)))
demography <- generate_pop_demography(N=N,times=times, age_min=180,prob_removal=1,aux=aux)

demography <- generate_pop_demography(N=N,times=times, age_min=180,prob_removal=1)

biomarker_map <- tibble(exposure_id=c("infection"),biomarker_id=c("IgG"))
reformat_biomarker_map(biomarker_map)

## 3. FOE models
foe_pars <- c(beta=0.1,gamma=1/7,I0=1/10000)

model_pars <- list(
    bind_rows(
        data.frame(i=1, t=10, x=1,biomarker_id=1,name=c("y0","y1","alpha","r","t1"), value=c(10,200, 0.07,2.16,5)),
        data.frame(i=1, t=60, x=1,biomarker_id=1,name=c("y1","alpha","r","t1"), value=c(250, 0.03,2.16,5))
    ))

y <- sapply(0:720, function(x) antibody_model_typhoid(1, x, 1, NULL, NULL, model_pars, NULL))


typhoid <- function(t, y0, y1, beta, r, t1){
    mu <- (1/t1) * log(y1/y0)
    y <- numeric(length(t))
    y[t<=t1] <- y0*exp(mu*t[t <= t1])
    alpha <- 1/(r-1)
    if(t <= t1){
        y <- y0*exp(mu*t)
    } else {
        tau <- t - t1
        y <- y1*(1+beta*tau)^(-alpha)
    }
    y
}

y1s <- rlnorm(N, log(200),1)
y0s <- rlnorm(N, log(13),0.75)
alphas <- rlnorm(N, log(0.07),0.75)
rs <- rlnorm(N, log(2.16), 0.1)
t1s <- rlnorm(N, log(5.5), 0.75)
times <- seq(0,365,by=1)
trajs <- matrix(NA, nrow=length(t), ncol=N)
for(i in 1:N){
    trajs[,i] <- sapply(times, function(t) typhoid(t, y0s[i],y1s[i],alphas[i],rs[i],t1s[i]))
}
mean_traj <- data.frame(t=t,y=apply(trajs, 1, median))
melted_trajs <- reshape2::melt(trajs)
colnames(melted_trajs) <- c("t","i","y")

ggplot(melted_trajs) + geom_line(aes(x=t,y=y,group=i),alpha=0.25) + 
    geom_line(data=mean_traj,aes(x=t,y=y),col="red") +
    scale_y_log10(limits=c(1,10000))

Rprof(tmp<-tempfile())
plot_antibody_model(antibody_model_monophasic, 5, model_pars=model_pars, biomarker_map=biomarker_map)
Rprof(NULL)
wow <- summaryRprof(tmp)
wow$by.self %>% arrange(-total.pct)



f <- function(){
    draw_parameters_fixed_fx(1,1,1,1,NULL,NULL,model_pars)
}
kinetics_pars <- list(draw_parameters_fixed_fx(1,1,1,1,NULL,NULL,model_pars) %>% drop_na())
exposure_history <- array(0,dim=c(10,50,2))
exposure_history[1,1,1] <- 1
f_ab <- function(){
    antibody_model_monophasic(1,25,1,exposure_history,NULL,kinetics_pars,biomarker_map)
}
f_ab()
microbenchmark::microbenchmark(f(), f_ab())

model_pars <- read.csv("~/Documents/GitHub/serosim/inst/extdata/model_pars_test_1.csv") %>% drop_na()
draw_parameters_random_fx(1,1,1,1,NULL,NULL,model_pars)
antibody_model_biphasic
y <- plot_antibody_model(antibody_model_biphasic, 5, model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx, biomarker_map=biomarker_map)


