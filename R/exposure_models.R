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
#' @description Probability of exposure depends on the force of exposure at the current time t for group g modulated by relevant demographic elements specified within the dem_mod. Within dem_mod, users can select which demographic elements affect the probability of exposure and by how much.
#' 
#' @param i Individual
#' @param t time
#' @param x exposure
#' @param g group
#' @param foe_pars A 3D array providing the force of exposure for each exposure ID, group and time.
#' @param demography A tibble of relevant demographic information for each individual in the simulation.
#' @param dem_mod A tibble specifying the modifier(how much each input affects probability of exposure) for each demographic elements; column names are column, value, modifier. Entries in column and value must match format in demography table. All column and value combinations in demography must have a modifier value within this tibble. Users can also add age modifier(how much each age affects probability of exposure). The column name will be "age" with the entry being individual's ages. 
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
  exps<-unique(dem_mod$exposure_id)
  
  ## If there are defined modifiers for this exposure then
  if(x  %in% exps){ 
    
    ## Pull out unique columns (demography variables) in dem_mod 
    cols<-unique(dem_mod$column)
    
    ## Pull individual's information within demography 
    individualnum<-i
    demography_tmp<-data.table(demography)
    demography_tmp<-demography_tmp[demography_tmp$i==individualnum & demography_tmp$times==t,]
    
    ## Convert dem_mod to data table
    mod2<-data.table(dem_mod)
    
    for (col in seq_along(cols)){ ## For each unique column entry in dem_mod 
      ## Pull the column name
      colname<-cols[col] 
      if(colname!="age"){
        ## Pull the individual's column entry within demography 
        entry<-demography_tmp[[colname]]
        ## Find the modifier within mod tibble 
        modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==entry]
        ## Multiply modifier by p 
        p <- p*modifier
      }
      if(colname=="age"){
        ## Calculate individual's current age
        ## Pull individual's birth time
        birth_time<-demography_tmp$birth
        curr_age<-floor((t-birth_time)/t_in_year)
        
        ## Find the age modifier within dem_mod tibble
        modifier<-mod2$modifier[mod2$exposure_id==x & mod2$column==colname & mod2$value==curr_age]
        ## Multiply modifier by p 
        p <- p*modifier
      }
    }
    p_exp<-1-exp(-p)
    p_exp
  }
}



  

