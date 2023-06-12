#list.of.packages <- c("ggplot2", "tidyverse","data.table","patchwork","reshape2", "tidyr","dplyr", "deSolve", "stats")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

#library(tidyverse)
#library(data.table)
#library(ggplot2)
#library(patchwork)
#library(reshape2)
#library(tidyr)
#library(dplyr)
#library(deSolve)
#library(stats)


## Define global variables needed to avoid "no visible binding for global variable" note 
utils::globalVariables(c('i','t','b','biomarker_id','titer','mean_titer','starting_biomarkers',
                          'boost_mod', 'biomarker', 'value', 'Time', 'Group', 
                         'Exposure', 'observed', 'times', 'x', 'birth', 'removal', 'group',
                         'exposure_id', "Exposure ID", 'cross_reactivity','curr_biom_quant', 'block'))
