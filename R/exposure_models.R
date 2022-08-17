#' Exposure Model Simple- Force of Infection
#' 
#' @description This is a simple exposure model where the probability of infection depends on the force of infection at that time for that group
#'
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param g group
#' @param lambdas Force of infection array 
#' @param demography Demography information 
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_simple_FOI <- function(i, t, e, g, lambdas, demography, ...){
  p <- lambdas[g, t, e]
  p_exp<-1-exp(-p)
  p_exp
}

#' Exposure Model Modified By Relevant Demographic Elements 
#'  
#' @description Probability of exposure depends on the force of infection at the current time t for group g modulated by relevant demographic elements specified within the mod. Within mod, users can select which demographic elements affect the probability of exposure and by how much.
#' 
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param g group
#' @param lambdas Force of infection array 
#' @param demography Demography information 
#' @param mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param ... 
#'  
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_dem_mod <- function(i, t, e, g, lambdas, demography, mod, ...){
  ## Find the force of infection 
  p <- lambdas[g, t, e]
  
  ## Pull out unique columns in mod 
  cols<-mod$column %>% unique()
  
  ## Pull individual's information within demography 
  individualnum<-i
  demography_tmp<-demography %>% filter(demography$i==individualnum & demography$times==t)
  
  for (col in seq_along(cols)){ ## For each unique column entry in mod 
    ## Pull the column name
    colname<-cols[col] 
    ## Pull the individual's column entry within demography 
    entry<-demography_tmp[colnames(demography_tmp)==colname]
    entry<-pull(entry)
    ## Find the modifier within mod tibble 
    modifier<-mod$modifier[mod$column==colname & mod$value==entry]
    ## Multiply modifier by p
    p <- p*modifier
  }
  p_exp<-1-exp(-p)
  p_exp
}


#' Exposure Model Modified By Age
#' 
#' @description Probability of exposure depends on the force of infection at the current time t for group g modulated by the individual's age specified within age_mod
#' 
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param g group
#' @param lambdas Force of infection array 
#' @param demography Demography information 
#' @param age_mod A tibble specifying the age modifier(how much each age affects probability of exposure); column names are age and modifier. All ages in simulation must have a modifier value within this tibble.
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_age_mod <- function(i, t, e, g, lambdas, demography, age_mod, ...){
  ## Find the force of infection 
  p <- lambdas[g, t, e]
  
  ## Calculate individual's current age
  curr_age<-t-birth_time ## Individual's birth time already gets pulled earlier on within the runserosim code so is it fine to just call it here?
  
  ## Find the modifier within age_mod tibble 
  age_modifier<-age_mod$modifier[age_mod$age==curr_age]
  
  ## Multiply modifier by p
  p <- p*age_modifier
  
  p_exp<-1-exp(-p)
  p_exp
}

#' Exposure Model Modified By Relevant Demographic Elements and Age
#' 
#' @description Probability of exposure depends on the force of infection at the current time t for group g modulated by relevant demographic elements (mod) and age (age_mod). This model combines exposure_model_V2 and exposure_model_V3
#' 
#' @param i Individual
#' @param t time
#' @param e exposure
#' @param g group
#' @param lambdas Force of infection array 
#' @param demography Demography information 
#' @param mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. 
#' @param age_mod A tibble specifying the age modifier(how much each age affects probability of exposure); column names are age and modifier. All ages in simulation must have a modifier value within this tibble.
#' @param ... 
#'
#' @return A probability of exposure is returned
#' @export
#'
#' @examples
exposure_model_dem_age_mod <- function(i, t, e, g, lambdas, demography, mod, age_mod, ...){
  ## Find the force of infection 
  p <- lambdas[g, t, e]
  
  ## Pull out unique columns in mod 
  cols<-mod$column %>% unique()
  
  ## Pull individual's information within demography 
  individualnum<-i
  demography_tmp<-demography %>% filter(demography$i==individualnum & demography$times==t)
  
  for (col in seq_along(cols)){ ## For each unique column entry in mod 
    ## Pull the column name
    colname<-cols[col] 
    ## Pull the individual's column entry within demography 
    entry<-demography_tmp[colnames(demography_tmp)==colname]
    entry<-pull(entry)
    ## Find the modifier within mod tibble 
    modifier<-mod$modifier[mod$column==colname & mod$value==entry]
    ## Multiply modifier by p
    p <- p*modifier
  }
  
  ## Calculate individual's current age
  curr_age<-t-birth_time ## Individual's birth time already gets pulled earlier on within the runserosim code so is it fine to just call it here?
  
  ## Find the age modifier within age_mod tibble 
  age_modifier<-age_mod$modifier[age_mod$age==curr_age]
  
  ## Multiply age modifier by p
  p <- p*age_modifier
  
  p_exp<-1-exp(-p)
  p_exp
}
  

