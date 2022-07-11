#' Observation Process
#'
#' @param titres The reshaped data set containing antibody titre for individuals at all time steps for each pathogen
#' @param obs_sd The distribution standard deviation; this represents the noise in the sampling process
#' @param distribution The distribution type; defaults to normal
#'
#' @return The titres data set is returned with a new column to represent the observed titre
#' @export
#'
#' @examples
observation_process<-function(titres, obs_sd, distribution="normal"){
  if(distribution=="normal"){
    titres$`Observed titre` <- stats::rnorm(nrow(titres),titres$Titre,obs_sd)
    return(titres)
  }
  if(distribution=="log-normal"){
    titres$`Observed titre` <- stats::rlnorm(nrow(titres),titres$Titre,obs_sd)
    return(titres)
  }
}
