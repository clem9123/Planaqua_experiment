
# importation des library
library(tidyverse)
library(readxl)

# importation des donn?es et mise en forme des colonnes
BDD_f <- read_excel("./data/data_final_ws_norm.xlsx", col_types = c("numeric","date","text","text","text","text","text","text","text","numeric","numeric","text","text"), na = "")
BDD_f <- BDD_f %>% mutate(Lake_capture = as.factor(Lake_capture),
                          Method_capture =as.factor(Method_capture),
                          Session_capture = as.factor(Session_capture),
                          Passage_net = as.factor(Passage_net),
                          Obs_status = as.factor(Obs_status),
                          Species = as.factor(Species),
                          Lake_released = as.factor(Lake_released))
summary(BDD_f)

# ordr des lacs en true false des traitements pour les graphes
BDD_f <- BDD_f %>% mutate (Lake_capture = factor(Lake_capture,levels = c(1,8,3,6,11,14,9,16,2,7,4,5,12,13,10,15)))
BDD_f <- BDD_f %>% mutate (Lake_released = factor(Lake_released,levels = c(1,8,3,6,11,14,9,16,2,7,4,5,12,13,10,15)))

# ajout de la colonne ann?e 
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

# cr?ation d'une colonne log_weight et Log_size
BDD_f <- BDD_f %>% mutate(Log_weight = log(Weight))
BDD_f <- BDD_f %>% mutate(Log_size = log(Size))

# Ajout de la colonne Tag_year
BDD_f <- BDD_f %>% group_by(Tag_id) %>% 
  mutate(Tag_year = ifelse(Tag_id %in%  c("juvenile","no_tag"), NA, min(as.character(Year)))) %>% 
  ungroup()

# creation du tableau roach adulte pris dans les trawls lors des sessions classiques
# ajout de la colonne annee tag


BDD_g <- BDD_f %>% filter(Species == "roach" & !(Method_capture %in% c("creel")) &
                            !(Session_capture %in% c("B","Z")))
BDD_g <- BDD_g %>% mutate (Obs_status = case_when(Tag_year %in% c("no_tag","juvenile")~ "capture",
                                                  Year == Tag_year ~ "capture",
                                                  Year != Tag_year ~ "recapture"))


BDD_p <- BDD_f %>% filter(Species == "perch" & !(Method_capture %in% c("creel")) &
                            !(Session_capture %in% c("B","Z")))
BDD_p <- BDD_p %>% mutate (Obs_status = case_when(Tag_year %in% c("no_tag","juvenile")~ "capture",
                                                  Year == Tag_year ~ "capture",
                                                  Year != Tag_year ~ "recapture"))

# cr?ation d'un tableau avec les nombres de poisson par ann?e et lac avec les poids et taille moyenne

BDD_n <- BDD_g %>% group_by(Year, Lake, Nutrients, Perch, Treatment) %>%
                     summarise(Mean_logweight = mean(Log_weight, na.rm =T),
                               Mean_logsize = mean (Log_size, na.rm =T),
                               Mean_weight = mean(Weight, na.rm =T),
                               Mean_size = mean (Size, na.rm =T),
                               nb_poisson = n(),
                               nb_juv = sum(Tag_id == "juvenile", na.rm = T),
                               biomasse = sum(Weight, na.rm =T)) %>%
                     ungroup()

BDD_a <- BDD_g %>% filter(Tag_id != "juvenile" & Tag_id != "no_tag")