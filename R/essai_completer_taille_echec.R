################## data

tr <- BDD_a %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
mooved <- tr[duplicated(tr$Tag_id),]$Tag_id
tr <- tr %>% filter ( !duplicated(tr$Tag_id))
tr <- tr %>% mutate (Treatment = ifelse(Tag_id %in% mooved, NA, Treatment))

la <- BDD_a %>% group_by(Tag_id, Lake) %>% summarize() %>% ungroup()
mooved <- la[duplicated(la$Tag_id),]$Tag_id
la <- la %>% filter ( !duplicated(la$Tag_id))
la <- la %>% mutate (Lake = ifelse(Tag_id %in% mooved, NA, Lake))

df <- BDD_a %>% ungroup () %>% pivot_wider(id_cols = Tag_id, names_from = Year, values_from = Size)
df <- df %>% merge(la) %>% merge(tr)

df <- df %>% filter (!is.na(Treatment)) %>% filter(!is.na(Lake))

set.seed(10)
df <- df[sample(1:nrow(df), 100),]

y <- df[2:7]
y[!is.na(y)] <- 1
y[is.na(y)] <- 0

z.inits <- function(ch){
  state <- ch
  state[state==0] <- 1
  #get.first <- function(x) min(which(x!=0))
  f <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    state[i,1:f[i]] <- NA
  }
  return(state)
}

get.first <- function(x) min(which(!is.na(x)))
f <- apply(y, 1, get.first)

get.first.size <- function(x) min(which(!is.na(x)))
size <- apply(df[2:7], 1, get.first.size)

jags.data <- list(y = y,
                  f = f,
                  Lake = df$Lake,
                  Treatment = df$Treatment,
                  Lake_treatment = Lake_treatment %>%
                    mutate(Lake = as.integer(Lake)) %>% 
                    mutate(Treatment = as.numeric(Treatment)) %>% select(Lake, Treatment),
                  nind = nrow(df),
                  noccas = 6,
                  ni = 1000,
                  size = size)

#### model de taille tout simple

model_taille_simple <- function(){
  # Likelihood
  for (i in 1:nind){
    for (t in (f[i]+1):noccas){
      y[i,t] ~ dnorm(y[i,t-1] + g*(270-y[i,t-1]), 1/(sd^2))
    }
  }
  # Prior
  g ~ dunif(0,100)
  sd ~ dunif(0,100)
  #Linf ~ dunif(100,300)
}

inits <- function(){
  list(p = runif(1,0,100), sd = runif(1,0,100))#, Linf =runif(1,100,300))
}

parameters = c("g","sd")#,"Linf")

Model_taille_simple <- jags.parallel(data = jags.data,
                                      inits = inits,
                                      parameters.to.save = parameters,
                                      model.file = model_taille_simple,
                                      n.chains = 4,
                                      n.iter = ni)

print(Model_taille_simple)
autocorr.plot(Model_taille_simple, ask = FALSE)
traceplot(Model_taille_simple, ask = F)

s_complet <- s
for (i in 1: nrow(s_complet)){
  for (t in (3: 8)){
    if (t != 8 & is.na(s_complet[i,t])){
      s_complet[i,t] <- s_complet[i,t-1] + 0.132*(s_complet[i,8]-s_complet[i,t-1])
    }
    if (t == 8 & is.na(s_complet[i,t])){
      s_complet[i,t] <- s_complet[i,t-1]
    }
  }
}
head(s_complet)

p_complet <- s_complet %>% pivot_longer(cols = 2:8)
colnames(p_complet) <- c("Tag_id", "Lake", "Treatment","Year", "Size")
head(p_complet)
ggplot(p_complet)+
  geom_line(aes(group = Tag_id, x = Year, y = Size, color = Treatment))



######### Model avec taille

cjs_taille <- function(){
  
  
  ####################
  ## Survie/Capture ##
  ####################
  
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t in (f[i]+1):noccas){
      z[i,t] ~ dbern(z[i,t-1]*phi)
      y[i,t] ~ dbern(z[i,t]*p)
    }
  }
  
  #Prior
  phi ~ dunif(0,1)
  p ~ dunif(0,1)
  
  ############
  ## Taille ##
  ############
  
  # Likelihood
  for (i in 1:nind){
    s[i,f[i]] <- size[i]
    for (t in (f[i]+1):noccas){
      s[i,t] ~ dnorm(s[i,t-1] + g, 1/(sd^2))
    }
  }
  
  # Prior
  g ~ dunif(0,100)
  sd ~ dunif(0,50)
}


inits <- list(list(phi = runif(1,0,1), p =  runif(1,0,1), g = runif(1,0,100), 
                   sd = runif(1,0,20), z = z.inits(y), s = df[2:7]),
              list(phi = runif(1,0,1), p =  runif(1,0,1), g = runif(1,0,100), 
                   sd = runif(1,0,20), z = z.inits(y), s = df[2:7]))


parameters = c("phi", "p", "g","sd")

CJS_taille <- jags(data = jags.data,
                   inits = inits,
                   parameters.to.save = parameters,
                   model.file = cjs_taille,
                   n.chains = 2,
                   n.iter = ni)







