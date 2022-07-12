#' Generate A Population Demography Data Set 
#'
#' @param N The number of individuals in the simulation
#' @param times The total number of time steps in the simulation
#' @param birth_times A vector of all individual's birth times
#' @param removal_times A vector of all individual's removal times
#' @param aux A list of the demography columns, the variable options and their distributions; defaults to NULL  
#'
#' @return
#' @export
#'
#' @examples
generate_pop_demography<-function(N, times, birth_times, removal_times, aux=NULL){
  if(is.null(aux)){
    df<- tibble(
      i=1:N,
      birth_times= birth_times,
      removal_times= removal_times)

    exp<- tidyr::expand_grid(1:N, times)
    exp<-dplyr::rename(exp,i="1:N")

    dem<- exp %>% dplyr::left_join(df, by="i")
    return(dem)
  }
  
  if(!is.null(aux)){
    vars <- NULL
    for(var in seq_along(aux)){
      vars[[var]] <- tibble(
        i=1:N,
        name=unlist(aux[[var]]["name"]),
        value=sample(aux[[var]][["options"]],size=N, prob=aux[[var]][["distribution"]],replace=TRUE)
      )
    }
    vars <- do.call("bind_rows",vars)
    vars <- vars %>% pivot_wider(names_from=name,values_from=value)
    
    df<- tibble(
      i=1:N,
      birth_times= birth_times,
      removal_times= removal_times)
    exp<- tidyr::expand_grid(1:N, times)
    exp<-dplyr::rename(exp,i="1:N")
    dem<- exp %>% dplyr::left_join(df, by="i")
    dem1<- dem %>% left_join(vars, by="i")
    return(dem1)
    
  }
}
  
  



  
  
  
  
  