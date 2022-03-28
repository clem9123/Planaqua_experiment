
# importation des library
library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(lme4)
library(GGally)
library(ggalluvial)
library(flexplot)
library(coefplot)
library(sjmisc)
library(sjPlot)
library(stats)

# importation des données et mise en forme des colonnes
setwd("C:/Users/sandrine/Documents/IEES_stage/code_et_donnee")
BDD_f <- read_excel("data/data_final_ws_norm.xlsx", col_types = c("numeric","date","text","text","text","text","text","text","text","numeric","numeric","text","text"), na = "")
BDD_f <- BDD_f %>% mutate(Lake_capture = as.factor(Lake_capture),
                          Method_capture =as.factor(Methode_capture),
                          Session_capture = as.factor(Session_capture),
                          Passage_net = as.factor(Passage_net),
                          Obs_status = as.factor(Obs_status),
                          Species = as.factor(Species),
                          Lake_released = as.factor(Lake_released))
summary(BDD_f)

# ordr des lacs en true false des traitements pour les graphes
BDD_f <- BDD_f %>% mutate (Lake_capture = factor(Lake_capture,levels = c(1,8,3,6,11,14,9,16,2,7,4,5,12,13,10,15)))
BDD_f <- BDD_f %>% mutate (Lake_released = factor(Lake_released,levels = c(1,8,3,6,11,14,9,16,2,7,4,5,12,13,10,15)))

# ajout de la colonne année 
BDD_f <- BDD_f %>% mutate (Year = as.factor(format(Date,"%Y")))

# ajout des colonnes de traitement des lacs : Nutrients, Perch, Treatment 
# en fonction du lac capture ou lac relacher si c'est une intro du lac relacher
Lake_treatment <- read_excel("data/Lake_treatment.xlsx", col_types =c("text","logical","logical","text"))
Lake_treatment <- Lake_treatment %>% mutate(Lake = as.factor(Lake),
                                            Treatment = as.factor(Treatment))
BDD_f <- BDD_f %>% mutate (Lake = case_when(is.na(Lake_capture) ~ Lake_released,
                                            TRUE ~ Lake_capture))
BDD_f <- merge (BDD_f, Lake_treatment)
BDD_f <- BDD_f %>% arrange(Index)

# création d'une colonne log_weight et Log_size
BDD_f <- BDD_f %>% mutate(Log_weight = log(Weight))
BDD_f <- BDD_f %>% mutate(Log_size = log(Size))

# création du tableau gardon adulte pris dans les chaluts lors des sessions classiques
# ajout de la colonne année tag
BDD_a <- BDD_f %>% filter(Species == "gardon" & 
                        Tag_id != "juvenile" &
                        !(Method_capture %in% c("Nasse")) &
                        !(Session_capture %in% c("B","Z")) &
                          (Year == "2016" | Obs_status != "introduction")) %>% 
                  group_by(Tag_id) %>% 
                  mutate(Tag_year = case_when(Tag_id == "juvenile" ~ "juvenile" ,
                                              is.na(Tag_id) ~ "NA",
                                              TRUE ~ min(as.character(Year))))
BDD_a$Tag_year <- BDD_a$Tag_year %>% na_if("NA")

## modification de la colonne Obs_status en ne prenant pas en compte les captures par Nasse
BDD_a <- BDD_a %>% mutate (Obs_status = case_when(Tag_year %in% c(NA,"juvenile")~ "capture",
                                                  Year == Tag_year ~ "capture",
                                                  Year != Tag_year ~ "recapture"))

# création du même tableau que ci dessus mais avec les juvéniles en plus
BDD_j <- BDD_f %>% filter(Species == "gardon" &
                          !(Method_capture %in% c("Nasse")) &
                          !(Session_capture %in% c("B","Z"))&
                          (Year == "2016" | Obs_status != "introduction")) %>% 
                    group_by(Tag_id) %>% 
                    mutate(Tag_year = case_when(Tag_id == "juvenile" ~ "juvenile" ,
                                      is.na(Tag_id) ~ "NA",
                                      TRUE ~ min(as.character(Year))))
BDD_j$Tag_year <- BDD_j$Tag_year %>% na_if("NA")

BDD_j <- BDD_j %>% mutate (Obs_status = case_when(Tag_year %in% c(NA,"juvenile") ~ "capture",
                                                  Year == Tag_year ~ "capture",
                                                  Year != Tag_year ~ "recapture"))

# création d'un tableau avec les nombres de poisson par année et lac avec les poids et taille moyenne

BDD_n <- BDD_a %>% group_by(Year, Lake, Nutrients, Perch, Treatment) %>%
                     summarise(Mean_logweight = mean(Log_weight, na.rm =T),
                               Mean_logsize = mean (Log_size, na.rm =T),
                               Mean_weight = mean(Weight, na.rm =T),
                               Mean_size = mean (Size, na.rm =T),
                               nb_poisson = n(),
                               nb_juv = sum(Tag_id == "juvenile", na.rm = T),
                               biomasse = sum(Weight, na.rm =T)) %>%
                     ungroup()




