list.of.packages <- c("ggplot2", "tidyverse","data.table","patchwork","reshape2", "tidyr","dplyr", "deSolve")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(reshape2)
library(tidyr)
library(dplyr)
library(deSolve)
