# TOUJOURS NECESSAIRE

library(R2jags)
library(tidyverse)
library(arules)
library(patchwork)

source("R/importation_data.R")

load("R/model_final/Model_Treatment.RData")
load("R/model_final/Model_Treatment_time_and_lake.RData")
load("R/model_final/Model_Treatment_time_x_lake.RData")
load("R/model_final/Model_Treatment_capture.RData")

file.edit("R/model_final/GRAPH_model_final.Rmd")

#######################

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

source("R/creation_tableau_taille_complet.R")

rm(tbl,t , t0, f, i, mooved, Linf, K, get.first, s, la, tr)