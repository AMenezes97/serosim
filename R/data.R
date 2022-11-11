#' Example demography tibble 
#'
#' Example of a demography tibble produced by generate_pop_demography and required as an input for runserosim
#' 
#' @docType data
#' @name example_demography
#' @usage data(example_demography)
#' @format A data frame with 12000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{times}{string numeric values for each time step in the simulation}
#'     \item{birth}{numeric value indicating that individual's birth time}
#'     \item{removal}{numeric value indicating that individual's removal time if any}
#' }
#' @family example_data
"example_demography"

#' Example biomarker map
#'
#' Example biomarker map. This biomarker map must be converted to the numeric version (example_biomarker_map_numeric) before being input into runserosim. This biomarker map is for a simulation with two different exposure events(infection and vaccination) both containing the same biomarker(IgG titer).
#' 
#' @docType data
#' @name example_biomarker_map
#' @usage data(example_biomarker_map)
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'     \item{exposure_ID}{name of each exposure type}
#'     \item{biomarker_ID}{name of each biomarkers present within each exposure type}
#' }
#' @family example_data
"example_biomarker_map"

#' Example numeric biomarker map
#'
#' Example numeric biomarker map required as an input for runserosim. This biomarker map is for a simulation with two different exposure events both containing the same biomarker.
#' 
#' @docType data
#' @name example_biomarker_map_numeric
#' @usage data(example_biomarker_map_numeric)
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'     \item{exposure_ID}{numeric values for each exposure type}
#'     \item{biomarker_ID}{numeric values for each biomarker present within each exposure type}
#' }
#' @family example_data
"example_biomarker_map_numeric"

#' Example force of exposure parameters (foe_pars) 
#'
#' Example force of exposure parameters required (foe_pars) as an input for runserosim. foe_pars is a 3D array providing the force of exposure for each exposure ID, group and time. This foe_pars argument is for a simulation with 120 time steps and two different exposure events (dimension 3)
#' 
#' @docType data
#' @name example_foe_pars
#' @usage data(example_foe_pars)
#' @format A 3 dimensional array containing the force of exposure at each time step for each exposure ID:
#' \describe{
#'     \item{dimension 1}{group dimension is 1  since no groups were specified in example_demography. Since there are no specified groups, all individuals are under the same force of exposure parameter}
#'     \item{dimenstion 2}{time step dimension is 120 since there are 120 time steps in the simulation}
#'     \item{dimenstion 3}{exposure ID dimension for 2 exposure types}
#' }
#' @family example_data
"example_foe_pars"

#' Example model parameters (model_pars)
#'
#' Example model parameters (model_pars). This model_pars input must be converted to the numeric version (example_model_pars_numeric) before being input into runserosim. This example model_pars is for a simulation two different exposure events both containing the same biomarker. model_pars argument is responsible for storing parameter information needed for the antibody model, observation model and immunity model.
#' In this example, model_pars has parameters for a monophasic antibody model and for an observational model with noise.
#'  
#' 
#' @docType data
#' @name example_model_pars
#' @usage data(example_model_pars)
#' @format A data frame with 5 rows and 6 variables:
#' \describe{
#'     \item{exposure_ID}{name of each exposure type present in example biomarker_map}
#'     \item{biomarker_ID}{name of each biomarker present within each exposure type present in example biomarker_map}
#'     \item{name}{names of model parameters}
#'     \item{mean}{numeric values for the true paramter means}
#'     \item{sd}{numeric values for the true paramter standard deviation}
#'     \item{distribution}{distribution type from which the draw_paramaters function will be sampling a parameter }
#' }
#' @family example_data
"example_model_pars"

#' Example numeric model parameters (model_pars)
#'
#' Example numeric model parameters (model_pars) required as an input for runserosim. This example model_pars is for a simulation two different exposure events both containing the same biomarker. model_pars argument is responsible for storing parameter information needed for the antibody model, observation model and immunity model.
#' In this example, model_pars has parameters for a monophasic antibody model and for an observational model with noise.
#'  
#' 
#' @docType data
#' @name example_model_pars_numeric
#' @usage data(example_model_pars_numeric)
#' @format A data frame with 5 rows and 6 variables:
#' \describe{
#'     \item{exposure_ID}{numeric values for each exposure type present in example biomarker_map}
#'     \item{biomarker_ID}{numeric values for each biomarker type present within each exposure type present in example biomarker_map}
#'     \item{name}{names of model parameters}
#'     \item{mean}{numeric values for the true paramter means}
#'     \item{sd}{numeric values for the true paramter standard deviation}
#'     \item{distribution}{distribution type from which the draw_paramaters function will be sampling a parameter }
#' }
#' @family example_data
"example_model_pars_numeric"

#' Example exposure histories 
#'
#' Example exposure history runserosim output. This example data set output contains exposure history for individuals at all time steps for each exposure ID
#' 
#' @docType data
#' @name example_exposure_histories
#' @usage data(example_exposure_histories)
#' @format A data frame with 24000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{e}{numeric values for each exposure type}
#'     \item{value}{binary value indicating if an exposure occured 1 or not 0}
#' }
#' @family example_data
"example_exposure_histories"

#' Example exposure probability 
#'
#' Example exposure probability runserosim output. This example data set output contains exposure probability for individuals at all time steps for each exposure ID
#' 
#' @docType data
#' @name example_exposure_probabilities
#' @usage data(example_exposure_probabilities)
#' @format A data frame with 24000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{e}{numeric values for each exposure type}
#'     \item{value}{numeric values indicating exposure probability}
#' }
#' @family example_data
"example_exposure_probabilities"

#' Example antibody states
#'
#' Example antibody states runserosim output. This example data set output contains antibody titers for individuals at all time steps for all biomarkers listed in example_biomarker_map
#' 
#' @docType data
#' @name example_antibody_states
#' @usage data(example_antibody_states)
#' @format A data frame with 12000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{b}{numeric values for each biomarker type; this example only has 1 biomarker type}
#'     \item{value}{numeric values indicating antibody titer}
#' }
#' @family example_data
"example_antibody_states"

#' Example observed antibody states
#'
#' Example observed antibody states runserosim output. This example data set output contains observed antibody titers for all individuals at the time of observation (t=120) for all biomarker 1
#' 
#' @docType data
#' @name example_observed_antibody_states
#' @usage data(example_observed_antibody_states)
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step observed}
#'     \item{b}{numeric values for each biomarker type; this example only has 1 biomarker type}
#'     \item{value}{numeric values indicating antibody titer}
#'     \item{observed}{numeric values indicating observed antibody titer}
#' }
#' @family example_data
"example_observed_antibody_states"