#' Example demography tibble 
#'
#' Example of a demography tibble produced by `generate_pop_demography` and required as an input for `runserosim`
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
#' Example biomarker map. This biomarker map must be converted to the numeric version (`example_biomarker_map_numeric`) before being input into `runserosim`. This biomarker map is for a simulation with two different exposure events (infection and vaccination) both containing the same biomarker (IgG titer).
#' 
#' @docType data
#' @name example_biomarker_map
#' @usage data(example_biomarker_map)
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'     \item{exposure_id}{name of each exposure type}
#'     \item{biomarker_id}{name of each biomarkers present within each exposure type}
#' }
#' @family example_data
"example_biomarker_map"

#' Example numeric biomarker map
#'
#' Example numeric biomarker map required as an input for `runserosim`. This biomarker map is for a simulation with two different exposure events both containing the same biomarker.
#' 
#' @docType data
#' @name example_biomarker_map_numeric
#' @usage data(example_biomarker_map_numeric)
#' @format A data frame with 2 rows and 2 variables:
#' \describe{
#'     \item{exposure_id}{numeric values for each exposure type}
#'     \item{biomarker_id}{numeric values for each biomarker present within each exposure type}
#' }
#' @family example_data
"example_biomarker_map_numeric"

#' Example force of exposure parameters (`foe_pars`) 
#'
#' Example force of exposure parameters required (`foe_pars`) as an input for `runserosim`. `foe_pars` is a 3D array providing the force of exposure for each exposure ID, group and time. This `foe_pars` argument is for a simulation with 120 time steps and two different exposure events (dimension 3)
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

#' Example model parameters (`model_pars`)
#'
#' Example model parameters (`model_pars`). This `model_pars` input must be converted to the numeric version (`example_model_pars_numeric`) before being input into `runserosim`. This example `model_pars` is for a simulation with two different exposure events both containing the same biomarker. `model_pars` argument is responsible for storing parameter information needed for the antibody model, observation model and immunity model.
#' In this example, `model_pars` has parameters for a monophasic antibody model and for an observational model with noise.
#'  
#' 
#' @docType data
#' @name example_model_pars
#' @usage data(example_model_pars)
#' @format A data frame with 5 rows and 6 variables:
#' \describe{
#'     \item{exposure_id}{name of each exposure type present in example biomarker_map}
#'     \item{biomarker_id}{name of each biomarker present within each exposure type present in example biomarker_map}
#'     \item{name}{names of model parameters}
#'     \item{mean}{numeric values for the true paramter means}
#'     \item{sd}{numeric values for the true paramter standard deviation}
#'     \item{distribution}{distribution type from which the draw_paramaters function will be sampling a parameter }
#' }
#' @family example_data
"example_model_pars"

#' Example model parameters for biphasic waning model
#'
#' This example `model_pars` is for a simulation three different exposure events (one vaccination, two infection types) corresponding to two biomarkers. `model_pars` argument is responsible for storing parameter information needed for the antibody model, observation model and immunity model. In this example, `model_pars` has parameters for a biphasic antibody model and for an observational model with noise. This `model_pars` input must be converted to the numeric version (`reformat_biomarker_map`) before being input into `runserosim`. 
#'  
#' 
#' @docType data
#' @name example_model_pars_biphasic
#' @usage data(example_model_pars_biphasic)
#' @format A data frame with 5 rows and 6 variables:
#' \describe{
#'     \item{exposure_id}{name of each exposure type present in example biomarker_map}
#'     \item{biomarker_id}{name of each biomarker present within each exposure type present in example biomarker_map}
#'     \item{name}{names of model parameters}
#'     \item{mean}{numeric values for the true paramter means}
#'     \item{sd}{numeric values for the true paramter standard deviation}
#'     \item{distribution}{distribution type from which the draw_paramaters function will be sampling a parameter }
#' }
#' @family example_data
"example_model_pars_biphasic"


#' Example numeric model parameters (`model_pars`)
#'
#' Example numeric model parameters (`model_pars`) required as an input for `runserosim`. This example `model_pars` is for a simulation two different exposure events both containing the same biomarker. `model_pars` argument is responsible for storing parameter information needed for the antibody model, observation model and immunity model.
#' In this example, `model_pars` has parameters for a monophasic antibody model and for an observational model with noise.
#'  
#' 
#' @docType data
#' @name example_model_pars_numeric
#' @usage data(example_model_pars_numeric)
#' @format A data frame with 5 rows and 6 variables:
#' \describe{
#'     \item{exposure_id}{numeric values for each exposure type present in example biomarker_map}
#'     \item{biomarker_id}{numeric values for each biomarker type present within each exposure type present in example biomarker_map}
#'     \item{name}{names of model parameters}
#'     \item{mean}{numeric values for the true paramter means}
#'     \item{sd}{numeric values for the true paramter standard deviation}
#'     \item{distribution}{distribution type from which the draw_paramaters function will be sampling a parameter }
#' }
#' @family example_data
"example_model_pars_numeric"

#' Example exposure histories 
#'
#' Example exposure history `runserosim` output. This example data set output contains exposure history for individuals at all time steps for each exposure ID
#' 
#' @docType data
#' @name example_exposure_histories
#' @usage data(example_exposure_histories)
#' @format A data frame with 24000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{x}{numeric values for each exposure type}
#'     \item{value}{binary value indicating if an exposure occured 1 or not 0}
#' }
#' @family example_data
"example_exposure_histories"


#' Example exposure histories (wide)
#'
#' Example exposure history `runserosim` output. This example data set output contains exposure history for individuals at all time steps for each exposure ID. This is identical to `example_exposure_histories`, but is converted into a 3D array giving the format expected by most antibody model functions.
#' 
#' @docType data
#' @name example_exposure_histories_wide
#' @usage data(example_exposure_histories_wide)
#' @format A 3D array with dimensions corresponding to 1) individual; 2) time; 3) exposure ID, and values giving the exposure state
#' @family example_data
"example_exposure_histories_wide"


#' Example exposure probability 
#'
#' Example exposure probability `runserosim` output. This example data set output contains exposure probability for individuals at all time steps for each exposure ID
#' 
#' @docType data
#' @name example_exposure_probabilities
#' @usage data(example_exposure_probabilities)
#' @format A data frame with 24000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{x}{numeric values for each exposure type}
#'     \item{value}{numeric values indicating exposure probability}
#' }
#' @family example_data
"example_exposure_probabilities"

#' Example force of exposure 
#'
#' Example force of exposure `runserosim` output. This example data set output contains the force of exposure for individuals at all time steps for each exposure ID
#' 
#' @docType data
#' @name example_exposure_force
#' @usage data(example_exposure_force)
#' @format A data frame with 24000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{x}{numeric values for each exposure type}
#'     \item{value}{numeric values indicating force of exposure}
#' }
#' @family example_data
"example_exposure_force"

#' Example biomarker states
#'
#' Example biomarker states from `runserosim` output. This example data set output contains biomarker quantities for individuals at all time steps for all biomarkers listed in `example_biomarker_map`
#' 
#' @docType data
#' @name example_biomarker_states
#' @usage data(example_biomarker_states)
#' @format A data frame with 12000 rows and 4 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step in the simulation}
#'     \item{b}{numeric values for each biomarker type; this example only has 1 biomarker type}
#'     \item{value}{numeric values indicating biomarker quantity}
#' }
#' @family example_data
"example_biomarker_states"


#' Example biomarker states (wide)
#'
#' Example biomarker states in wide format from `runserosim` output. This is identical to `example_biomarker_states`, but is converted into a 3D array giving the format expected by most antibody model functions. This example data set output contains biomarker quantities for individuals at all time steps for all biomarkers listed in `example_biomarker_map`
#' 
#' @docType data
#' @name example_biomarker_states_wide
#' @usage data(example_biomarker_states_wide)
#' @format A 3D array with dimensions corresponding to 1) individual; 2) time; 3) biomarker, and values giving the latent biomarker quantity
#' @family example_data
"example_biomarker_states_wide"

#' Example observed biomarker states
#'
#' Example observed biomarker states `runserosim` output. This example data set output contains observed biomarker quantities for all individuals at the time of observation (t=120) for all biomarker 1
#' 
#' @docType data
#' @name example_observed_biomarker_states
#' @usage data(example_observed_biomarker_states)
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'     \item{i}{string numeric values for each individual}
#'     \item{t}{string numeric values for each time step observed}
#'     \item{b}{numeric values for each biomarker type; this example only has 1 biomarker type}
#'     \item{value}{numeric values indicating biomarker quantity}
#'     \item{observed}{numeric values indicating observed biomarker quantity}
#' }
#' @family example_data
"example_observed_biomarker_states"