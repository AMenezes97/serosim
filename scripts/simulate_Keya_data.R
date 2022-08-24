## Simulate a data set to input to Keya's latent class model 

## Data structure
## mom_indicator:  binary indicator of if an individual had mother recall of any mcv dose among children with or without card
## card_indicator: binary indicator of if an individual had record (via card) of any mcv dose
## none_indicator = binary indicator of if an individual had neither card nor recall
## boost_indicator = binary indicator of if an individual has a suspected boost
## titer = log antibody concentrations of each individual 
## age = age in years of each individual. Age here is a double (e.g. 4.39) (Age ranges 0.7 to 14.99)
## ageyr = age in years of each individual as an integer (floor of age) 
## ageyr_obs = age in integers among those with observed vaccine information
## ageyr_mis = age in integers among those with missing vaccine information 
## location = one of 4 pre-campaign locations an individual is from 


## Figure out which, if any,  demography elements will get added prior to the simulation 
## and what their affect will be on any of the subsequent models

## Figure out which variables will get assigned post simulation (mother recall/card indicator)
## These variables will depend on an individual's vaccination and exposure history

## For example: If vacc==1 then some proportion of those individuals have have a 1 for card indicator and mother recall
## We then feed this information into the model and have it recall these numbers

## James said we should simulate the full data before Keya manipulated and removed certain variables. 
## Does this mean we simulate the data set before momrec and docrec were combined into the overall recall variables?
