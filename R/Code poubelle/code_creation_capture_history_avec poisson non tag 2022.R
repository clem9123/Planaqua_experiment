subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Lake
data.frame(Tag_id = "no_tag", 
           X2016 = 0, X2017 = 0, X2018 = 0, X2019 = 0,X2020 = 0, X2021 = 0, X2022 = 1, 
           Treatment = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Treatment, 
           Lake = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Lake )
data.frame(Tag_id = "no_tag", X2016 = 0, X2017 = 0, X2018 = 0, X2019 = 0,X2020 = 0, X2021 = 0,
           X2022 = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Size, 
           Treatment = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Treatment,
           Lake = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Lake )

################# data avec 2022 et les roachs no_tag de 2022
s <- s2

tr <- BDD_a %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

t <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from =Lake)
t <- merge(t,tr)
t <- merge(t,la)

t <- data.frame(t$Tag_id, apply(t[2:8],2, function(x) ifelse(is.na(x),0,1)), t$Treatment, t$Lake)

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

colnames(t)<- c("Tag_id", "X2016","X2017","X2018","X2019","X2020","X2021","X2022", "Treatment", "Lake")
t <- t %>% rbind(data.frame(Tag_id = "no_tag", 
                            X2016 = 0, X2017 = 0, X2018 = 0, X2019 = 0,X2020 = 0, X2021 = 0, X2022 = 1, 
                            Treatment = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Treatment, 
                            Lake = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Lake ))
colnames(s)<- c("Tag_id", "X2016","X2017","X2018","X2019","X2020","X2021","X2022", "Treatment", "Lake")
s <- s %>% rbind(data.frame(Tag_id = "no_tag", X2016 = 0, X2017 = 0, X2018 = 0, X2019 = 0,X2020 = 0, X2021 = 0,
                            X2022 = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Size, 
                            Treatment = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Treatment,
                            Lake = subset(BDD_g, Tag_id=="no_tag" & Weight > 8)$Lake ))


#t <- na.omit(t) necessaire si j'utilise les lac ou  lest traitements

#set.seed(10)
#t <- t[sample(1:nrow(t), 100),]


CH <- as.matrix(t[2:8])

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

z.inits <- function(ch){
  state <- ch
  state[state==0] <- 1
  get.first <- function(x) min(which(x!=0))
  f <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    state[i,1:(f[i])] <- NA
  }
  return(state)
}

jags.data <- list(y = CH,
                  size = as.matrix(s[2:7]),
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 500,
                  zi = z.inits(CH))
