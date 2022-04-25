# necessite juste BDD_a

library(R2jags)

load("R/object/Model_g_Kla_Lla_2016.RData")

# remplissage du tableau des tailles
# création d'un tableau d'augmentation de taille par lac et année
tbl = data.frame()
t0 <- Model_g_Kla_Lla_2016[[2]]$mean$t0
for (i in 1:16){
  K <- Model_g_Kla_Lla_2016[[2]]$mean$K[i]
  Linf <- Model_g_Kla_Lla_2016[[2]]$mean$Linf[i]
  
  for (t in 2017:2021){
    tbl <- rbind(tbl, data.frame(i,t,Linf*(1-exp(-K*(t-t0)))-Linf*(1-exp(-K*(t-1-t0)))))
  }
}
colnames(tbl) <- c("Lake", "Year", "Size")
tbl <- tbl %>% mutate(Lake = factor(Lake), Year = factor(Year))
head(tbl)

ggplot(tbl)+
  geom_point(aes(x = Year, y = Size))+
  facet_wrap(~ Lake)

# Maintenant remplissage
tr <- BDD_a %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

s <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)
s <- s %>% merge(la) %>% merge(tr)



head(s)

get.first <- function(x) min(which(!is.na(x)))
f <- apply(s[2:8], 1, get.first)

for (i in 1:nrow(s)){
  if (!is.na(s$Lake[i])){
    t0 <- f[i]+1
    if (t0 <7){
      for (t in (t0+1):7){ # ou 7 je ne sais pas encore comment je preends en comte 2022
        if (is.na(s[i,t])){
          s[i,t] <- s[i,t-1] + subset(tbl, Lake == s$Lake[i] & Year == 2015+t-f[i])$Size
        }
      }
    }
    if (is.na(s[i,8])){s[i,8] <- s[i,7]}
  }
}

