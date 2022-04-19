# Tu peux récupérér BDD_a dans R/Object
# Le but de ce plot était de voir les mouvements des poissons en 2018

library(ggalluvial) # et cette page m'a pas mal aidé aussi : http://corybrunson.github.io/ggalluvial/articles/ggalluvial.html

# Vision des changements de lac en 2018 (sans les introductions) en alluvial plot


mooved_2018 <- BDD_p %>% filter(Year == "2018") %>% filter (!(Tag_id %in% c("no_tag","juvenile")))
mooved_2018 <- mooved_2018 %>% mutate(mooved = case_when(
  Lake_capture == Lake_released ~FALSE,
  TRUE ~TRUE))

mooved_2018 <- mooved_2018 %>% mutate(Lake_capture = factor(Lake_capture, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15))) %>%
                               mutate(Lake_released = factor(Lake_released, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15)))


ggplot(mooved_2018,
       aes(axis1 = Lake_capture, axis2 = Lake_released)) +
  geom_alluvium(aes(fill = mooved),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Lake1", "Lake2"), expand = c(.05, .05))


# pour 2019

mooved_2019 <- BDD_p %>% filter(Year == "2019") %>% filter (!(Tag_id %in% c("no_tag","juvenile")))
mooved_2019 <- mooved_2019 %>% mutate(mooved = case_when(
  Lake_capture == Lake_released ~FALSE,
  TRUE ~TRUE))

mooved_2019 <- mooved_2019 %>% mutate(Lake_capture = factor(Lake_capture, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15))) %>%
  mutate(Lake_released = factor(Lake_released, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15)))


ggplot(mooved_2019,
       aes(axis1 = Lake_capture, axis2 = Lake_released)) +
  geom_alluvium(aes(fill = mooved),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Lake1", "Lake2"), expand = c(.05, .05))

# pour 2020

mooved_2020 <- BDD_p %>% filter(Year == "2020") %>% filter (!(Tag_id %in% c("no_tag","juvenile")))
mooved_2020 <- mooved_2020 %>% mutate(mooved = case_when(
  Lake_capture == Lake_released ~FALSE,
  TRUE ~TRUE))

mooved_2020 <- mooved_2020 %>% mutate(Lake_capture = factor(Lake_capture, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15))) %>%
  mutate(Lake_released = factor(Lake_released, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15)))


ggplot(mooved_2020,
       aes(axis1 = Lake_capture, axis2 = Lake_released)) +
  geom_alluvium(aes(fill = mooved),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Lake1", "Lake2"), expand = c(.05, .05))

# pour 2021

mooved_2021 <- BDD_p %>% filter(Year == "2021") %>% filter (!(Tag_id %in% c("no_tag","juvenile")))
mooved_2021 <- mooved_2021 %>% mutate(mooved = case_when(
  Lake_capture == Lake_released ~FALSE,
  TRUE ~TRUE))

mooved_2021 <- mooved_2021 %>% mutate(Lake_capture = factor(Lake_capture, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15))) %>%
  mutate(Lake_released = factor(Lake_released, levels = c(1,8,11,14,2,7,12,13,3,6,9,16,4,5,10,15)))


ggplot(mooved_2021,
       aes(axis1 = Lake_capture, axis2 = Lake_released)) +
  geom_alluvium(aes(fill = mooved),width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum)))+
  scale_x_discrete(limits = c("Lake1", "Lake2"), expand = c(.05, .05))

