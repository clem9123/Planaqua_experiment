# Tu peux récupérér BDD_a dans R/Object
# Le but de ce plot était de voir les mouvements des poissons en 2018

library(ggplot2)
library(dplyr)
library(ggalluvial) # et cette page m'a pas mal aidé aussi : http://corybrunson.github.io/ggalluvial/articles/ggalluvial.html

# Vision des changements de lac en 2018 (sans les introductions) en alluvial plot

mooved_2018 <- BDD_a %>% filter(Year == "2018" & Obs_status != "introduction")
mooved_2018 <- mooved_2018 %>% mutate(mooved = case_when(
  Lake_capture == Lake_released ~FALSE,
  Lake_capture != Lake_released ~TRUE))


ggplot(mooved_2018,
       aes(axis1 = Lake_capture, axis2 = Lake_released)) +
  geom_alluvium(aes(fill = mooved),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "grey") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Lake1", "Lake2"), expand = c(.05, .05))
