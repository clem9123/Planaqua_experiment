library(R2jags)
library(tidyverse)
library(arules)

source("R/importation_data.R")
s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

source("R/creation_tableau_taille_complet.R")

rm(tbl,t , t0, f, i, mooved, Linf, K, get.first, s, la, tr)
