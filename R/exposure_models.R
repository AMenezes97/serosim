#' Exposure Model Simple- Force of Exposure (FOE)
#' 
#' @description This is a simple exposure model where the probability of exposure depends on the force of exposure at that time for that group
#'
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_simple_FOE <- function(i, t, x, g, foe_pars, demography, ...){
  p <- foe_pars[g, t, x]
  p_exp<-1-exp(-p)
  p_exp
}

#' Exposure Model Modified By Relevant Demographic Elements 
#'  
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements specified within the mod. Within mod, users can select which demographic elements affect the probability of exposure and by how much.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param ... 
#'  
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_dem_mod <- function(i, t, x, g, foe_pars, demography, dem_mod, ...){
  ## Find the force of exposure
  p <- foe_pars[g, t, x]
  
  ## Find which exposure events have modifiers
  exps<-dem_mod$exposure_id %>% unique()
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
  
  ## Pull out unique columns (demography variables) in dem_mod 
  cols<-dem_mod$column %>% unique()
  
  ## Pull individual's information within demography 
  individualnum<-i
  demography_tmp<-data.table(demography)
  demography_tmp<-demography_tmp[demography_tmp$i==individualnum & times==t,]
  
  for (col in seq_along(cols)){ ## For each unique column entry in dem_mod 
    ## Pull the column name
    colname<-cols[col] 
    ## Pull the individual's column entry within demography 
    entry<-demography_tmp[[colname]]
    ## Find the modifier within mod tibble 
    mod2<-data.table(dem_mod)
    modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
    ## Multiply modifier by p 
    p <- p*modifier
  }
  }
  p_exp<-1-exp(-p)
  p_exp
}


#' Exposure Model Modified By Age
#' 
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by the individual's age specified within age_mod
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param age_mod A tibble specifying the age modifier(how much each age affects probability of exposure); column names are age and modifier. All ages in simulation must have a modifier value within this tibble.
#' @param t_in_year The number of time steps in a year; defaults to 1 
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_age_mod <- function(i, t, x, g, foe_pars, demography, age_mod, t_in_year=1, ...){
  ## Find the force of exposure
  p <- foe_pars[g, t, x]
  
  ## Check if age modifiers apply to this particular exposure 
  exps_age<-age_mod$exposure_id %>% unique()
  
  ## If there are defined modifiers for that exposures then
  if(x  %in% exps_age){ 
  
    age_mod2<-data.table(age_mod)
    age_mod_tmp<-age_mod2[age_mod$exposure_id==x,]
    
  ## Calculate individual's current age
    ## Pull individual's information within demography 
    individualnum<-i
    demography_tmp<-data.table(demography)
    demography_tmp<-demography_tmp[demography_tmp$i==individualnum & times==t,]
    birth_time<-demography_tmp$birth
  curr_age<-floor((t-birth_time)/t_in_year) ## Individual's birth time already gets pulled earlier on within the runserosim code so is it fine to just call it here?
  
  ## Find the modifier within age_mod tibble 
  age_modifier<-age_mod_tmp$modifier[age_mod_tmp$age==curr_age]
  
  ## Multiply modifier by p
  p <- p*age_modifier
  }
  p_exp<-1-exp(-p)
  p_exp

}

#' Exposure Model Modified By Relevant Demographic Elements and Age
#' 
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements (mod) and age (age_mod). This model combines exposure_model_V2 and exposure_model_V3
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param age_mod A tibble specifying the age modifier(how much each age affects probability of exposure); column names are age and modifier. All ages in simulation must have a modifier value within this tibble.
#' @param t_in_year The number of time steps in a year; defaults to 1 
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_dem_age_mod <- function(i, t, x, g, foe_pars, demography, dem_mod, age_mod, t_in_year=1, ...){
  ## Find the force of exposure 
  p <- foe_pars[g, t, x]
  
  ## Find which exposure events have modifiers
  exps<-dem_mod$exposure_id %>% unique()
  
  ## Pull individual's information within demography 
  individualnum<-i
  demography_tmp<-data.table(demography)
  demography_tmp<-demography_tmp[demography_tmp$i==individualnum & times==t,]
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
    
    ## Pull out unique columns in den_mod 
    cols<-dem_mod$column %>% unique()
    
    for (column in seq_along(cols)){ ## For each unique column entry in dem_mod 
      ## Pull the column name
      colname<-cols[column] 
      ## Pull the individual's column entry within demography 
      entry<-demography_tmp[[colname]]
      ## Find the modifier within mod tibble 
      mod2<-data.table(dem_mod)
      modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
      ## Multiply modifier by p
      p <- p*modifier
    }
  }
  
  ## Check if age modifiers apply to this particular exposure 
  exps_age<-age_mod$exposure_id %>% unique()
  
  ## If there are defined modifiers for that exposures then
  if(x  %in% exps_age){ 
    
    age_mod2<-data.table(age_mod)
    age_mod_tmp<-age_mod2[age_mod$exposure_id==x,]
    
    ## Calculate individual's current age
    birth_time<-demography_tmp$birth
    curr_age<-floor((t-birth_time)/t_in_year) ## Individual's birth time already gets pulled earlier on within the runserosim code so is it fine to just call it here?
    
    ## Find the age modifier within age_mod tibble 
    age_modifier<-age_mod_tmp$modifier[age_mod_tmp$age==curr_age]
    
    ## Multiply age modifier by p
    p <- p*age_modifier
  }
    p_exp<-1-exp(-p)
    p_exp
  
}
  

