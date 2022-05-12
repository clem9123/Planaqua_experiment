# TOUJOURS NECESSAIRE

library(R2jags)
library(tidyverse)
library(arules)
library(patchwork)

source("R/importation_data.R")

load("R/object/Model_multi_treatment_time_abundance.RData")
load("R/object/Model_multi_treatment_corrected_abundance.RData")
load("R/object/Model_multi_treatment_random_corrected_abundance.RData")


load("R/object/Model_lake_time_abundance.RData") #Model_multi_lake_time_abundance
load("R/object/Model_multi_lake_corrected_abundance.RData")


#######################

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

source("R/creation_tableau_taille_complet.R")

rm(tbl,t , t0, f, i, mooved, Linf, K, get.first, s, la, tr)