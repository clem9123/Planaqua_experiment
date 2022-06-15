# -----------------
# Necessary library
# -----------------

library(R2jags)
library(tidyverse)
library(arules)
library(patchwork)
library(ggalluvial)
library(LaplacesDemon)
library(lme4)
library(scales)

# ----------------
# Data importation
# ----------------

source("R/importation_data.R")

#------------------------
# JAGS data for the model
#------------------------

# capture history data frame with size
s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)

# extraction of treatment covariate
tr <- BDD_a %>% filter(Year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup() #%>% filter(Year != "2022") 
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

# extraction of lake covariate
la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

# merge capture history and covariate
s <- merge(s,tr)
s <- merge(s,la)

# elimination of roaches who changed lake
s <- s %>% filter(!is.na(Lake))

# discretization of size in 3 size class
size_break <- c(0,120,180,300)
s_group <- data.frame(apply(s[2:7], 2, 
                            function(x) discretize(x,method = "fixed",
                                                   breaks = size_break, 
                                                   labels= c(1:(length(size_break)-1)))))
s_group [] <- apply(s_group, 2, as.numeric)

# extraction of size class at first capture
fs <- apply(s_group, 1, function(x) first(na.omit(x)))

s_group <- data.frame(apply(s_group,2, function(x) ifelse(is.na(x),4,x)))

# extraction of time of first capture
get.first <- function(x) min(which(x!=4))
f <- apply(s_group, 1, get.first)

# creation zinits (real state initial)
zinit <- s_group
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {
    if (j > f[i] & s_group[i,j]==4) {zinit[i,j] <- ifelse(zinit[i,j-1]==1,2,zinit[i,j-1])}
  }
}
for (i in 1:nrow(s_group)) {
  for (j in 1:ncol(s_group)) {    
    if (j <= f[i]) {zinit[i,j] <- NA}
  }
}
for (i in 1:nrow(s_group)) {
  if (f[i]<5){
    for (j in (f[i]+2):ncol(zinit)) {
      if (zinit[i,j]<zinit[i,j-1]) {zinit[i,j] <- zinit[i,j-1]}
    }
  }}
zinit <- as.matrix(zinit)

# data for the model
jags.data <- list(y = s_group,
                  zi = zinit,
                  f = f,
                  fs = fs,
                  nind = nrow(s_group),
                  noccas = ncol(s_group),
                  ni = 5000,
                  Lake = s$Lake,
                  Treatment = s$Treatment,
                  Tr1 = which(s$Treatment == 1),
                  Tr2 = which(s$Treatment == 2), 
                  Tr3 = which(s$Treatment == 3), 
                  Tr4 = which(s$Treatment == 4),
                  LT = Lake_treatment$Treatment)

