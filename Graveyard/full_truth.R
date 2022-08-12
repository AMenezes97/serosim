#' Create Serosurvey Observation Data Set With Individual's Infection and Vaccination History
#'
#' @param titres Titre data set
#' @param obs_time The observation time when titres were measured
#' @param infection_histories Infection history data set
#' @param vaccine_histories_reshaped Vaccine history data set
#' @param age Individual's age data set
#' @param N_pathogens The number of pathogens in the simulation
#'
#' @return A new data set is returned with individual's true and observed titres at the time of the serosurvey along with whether they were infected, vaccinated or both
#' @export
#'
#' @examples
full_truth <-function(titres, obs_time, infection_histories, vaccine_histories_reshaped, age, N_pathogens){
  if(N_pathogens==2){
    overallinfectedP1 <- infection_histories %>%
      dplyr::filter(Time<obs_time) %>%  #individuals who are infected at the last time step show up with a 0 titre
      dplyr::filter(Pathogen==1) %>%
      dplyr::select(Individual, 'Infected?', Pathogen) %>%
      dplyr::filter(`Infected?`==1)

    titres1 <- titres %>%
      dplyr::select(Individual,Time,Titre, Pathogen,`Observed titre`) %>%
      dplyr::left_join(overallinfectedP1, by="Individual")

    overallinfectedP2 <- infection_histories %>%
      dplyr::filter(Time<obs_time) %>%  #individuals who are infected at the last time step show up with a 0 titre
      dplyr::filter(Pathogen==2) %>%
      dplyr::select(Individual, 'Infected?', Pathogen) %>%
      dplyr::filter(`Infected?`==1)

    titres2 <- titres1 %>%
      dplyr::select(Individual,Time,Titre, Pathogen.x,`Observed titre`,`Infected?`,Pathogen.y) %>%
      dplyr::left_join(overallinfectedP2, by="Individual")

    overallvaccinated<- vaccine_histories_reshaped  %>%
      dplyr::filter(Time<obs_time) %>%
      dplyr::select(Individual, 'Vaccinated?') %>%
      dplyr::filter(`Vaccinated?`==1)

    titres3 <- titres2 %>%
      dplyr::select(Individual,Time,Titre,Pathogen.x,`Observed titre`,`Infected?.x`,Pathogen.y,`Infected?.y`, Pathogen) %>%
      dplyr::left_join(overallvaccinated, by="Individual")

    titres3$Infected1<-as.numeric(titres3$`Infected?.x`)
    titres3$Infected2<-as.numeric(titres3$`Infected?.y`)
    titres3$Vaccinated<-as.numeric(titres3$`Vaccinated?`)

    Events <- rep("Neither", times=nrow(titres3))
    Events[is.na(titres3$Infected1) & is.na(titres3$Infected2) & is.na(titres3$Vaccinated)] <- "No infections or vaccination"
    Events[titres3$Infected1==2 & is.na(titres3$Infected2) & is.na(titres3$Vaccinated)] <- "Infected1"
    Events[is.na(titres3$Infected1) & titres3$Infected2==2 & is.na(titres3$Vaccinated)] <- "Infected2"
    Events[titres3$Infected1==2 & titres3$Infected2==2  & is.na(titres3$Vaccinated)] <- "Infected12"
    Events[is.na(titres3$Infected1) & is.na(titres3$Infected2) & titres3$Vaccinated==2] <- "Vaccinated"
    Events[titres3$Infected1==2 & is.na(titres3$Infected2) & titres3$Vaccinated==2] <- "Infected1/Vacc"
    Events[is.na(titres3$Infected1) & titres3$Infected2==2 & titres3$Vaccinated==2] <- "Infected2/Vacc"
    Events[titres3$Infected1==2 & titres3$Infected2==2 & titres3$Vaccinated==2] <- "Infected12/Vacc"
    titres3$Events <- Events

    titres3$`Infected?.x` <-NULL
    titres3$Pathogen.y  <-NULL
    titres3$`Infected?.y` <-NULL
    titres3$Pathogen <-NULL
    titres3$`Vaccinated?` <-NULL

    serosurveyage<- age %>%
      dplyr::filter(Time==obs_time) %>%
      dplyr::filter(Pathogen==1) %>%
      dplyr::select(Individual, Age)


    serosurveyage$Age <-as.numeric(serosurveyage$Age)-1

    #Convert age to years
    serosurveyage <- serosurveyage %>%
      dplyr::mutate(Age=Age/12)

    titresobs <- titres3 %>%
      dplyr::filter(Time==obs_time) %>%
      dplyr::select(Individual,Time,Titre,Pathogen.x,`Observed titre`, Infected1, Infected2, Vaccinated, Events) %>%
      dplyr::left_join(serosurveyage, by="Individual")

    #Add age categories to titresobs
    AgeCat <- rep(NA,nrow(titresobs))
    AgeCat[titresobs$Age <1] <- "0"
    AgeCat[titresobs$Age >=1 & titresobs$Age < 2] <- "1"
    AgeCat[titresobs$Age >=2 & titresobs$Age < 3] <- "2"
    AgeCat[titresobs$Age >=3 & titresobs$Age < 4] <- "3"
    AgeCat[titresobs$Age >=4 & titresobs$Age < 5] <- "4"
    AgeCat[titresobs$Age >=5 & titresobs$Age < 6] <- "5"
    AgeCat[titresobs$Age >=6 & titresobs$Age < 7] <- "6"
    AgeCat[titresobs$Age >=7 & titresobs$Age < 8] <- "7"
    AgeCat[titresobs$Age >=8 & titresobs$Age < 9] <- "8"
    AgeCat[titresobs$Age >=9 & titresobs$Age < 10] <- "9"
    titresobs$AgeCat <- as.ordered(AgeCat)

    titresobs<-rename(titresobs, Pathogen=Pathogen.x)

    return(titresobs)
  }
  if(N_pathogens==1){

    overallinfected <- infection_histories %>%
      filter(Time<obs_time) %>% #individuals who are infected at the last time step show up with a 0 titre
      select(Individual, 'Infected?') %>%
      filter(`Infected?`==1)

    titres1 <- titres %>%
      select(Individual,Time,Titre, Pathogen,`Observed titre`) %>%
      left_join(overallinfected, by="Individual")

    overallvaccinated<- vaccine_histories_reshaped %>%
      filter(Time<obs_time) %>%
      select(Individual, 'Vaccinated?') %>%
      filter(`Vaccinated?`==1)

    titres2 <- titres1 %>%
      select(Individual,Time,Titre,Pathogen,`Observed titre`, 'Infected?') %>%
      left_join(overallvaccinated, by="Individual")

    titres2$InfectedNumeric<-as.numeric(titres2$`Infected?`)
    titres2$VaccinatedNumeric<-as.numeric(titres2$`Vaccinated?`)

    Events <- rep("Neither", times=nrow(titres2))
    Events[titres2$InfectedNumeric==2 & is.na(titres2$VaccinatedNumeric)] <- "Infected"
    Events[titres2$InfectedNumeric==2 & titres2$VaccinatedNumeric==2] <- "Both"
    Events[is.na(titres2$InfectedNumeric) & titres2$VaccinatedNumeric==2] <- "Vaccinated"
    titres2$Events <- Events

    titres2$`Infected?`<-NULL
    titres2$`Vaccinated?`<-NULL

    serosurveyage<- age %>%
      filter(Time==obs_time) %>%
      select(Individual, Age)

    serosurveyage$Age <-as.numeric(serosurveyage$Age)-1

    #Convert age to years
    serosurveyage <- serosurveyage %>%
      mutate(Age=Age/12)

    titresobs <- titres2 %>%
      filter(Time==obs_time) %>%
      select(Individual,Time,Titre,Pathogen,`Observed titre`, InfectedNumeric, VaccinatedNumeric, Events) %>%
      left_join(serosurveyage, by="Individual")


    #Add age categories to titresobs
    AgeCat <- rep(NA,nrow(titresobs))
    AgeCat[titresobs$Age <1] <- "0"
    AgeCat[titresobs$Age >=1 & titresobs$Age < 2] <- "1"
    AgeCat[titresobs$Age >=2 & titresobs$Age < 3] <- "2"
    AgeCat[titresobs$Age >=3 & titresobs$Age < 4] <- "3"
    AgeCat[titresobs$Age >=4 & titresobs$Age < 5] <- "4"
    AgeCat[titresobs$Age >=5 & titresobs$Age < 6] <- "5"
    AgeCat[titresobs$Age >=6 & titresobs$Age < 7] <- "6"
    AgeCat[titresobs$Age >=7 & titresobs$Age < 8] <- "7"
    AgeCat[titresobs$Age >=8 & titresobs$Age < 9] <- "8"
    AgeCat[titresobs$Age >=9 & titresobs$Age < 10] <- "9"
    titresobs$AgeCat <- as.ordered(AgeCat)
    return(titresobs)
  }

  }



