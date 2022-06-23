################# data
s <- taillle

tr <- BDD_a %>% filter (Year != 2022) %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
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
t <- na.omit(t)

set.seed(11)
t <- t[sample(1:nrow(t), 500),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

CH <- as.matrix(t[2:7])

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

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

jags.data <- list(y = CH,
                  size = as.matrix(s[2:7]),
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 1000,
                  zi = z.inits(CH))

##################### Model 1
# avec b qui pourra varier avec les traitements pour phi

# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a = runif(1,-0.2,0), b = runif(1,-0.2,0.2), z = zi)
}

parameters = c("a","b")

cjs_taille <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a*size[i,t]
      logit(phi[i,t]) <- b*size[i,t]
      z[i,t] ~ dbern(p[i,t]*z[i,t-1])
      y[i,t] ~ dbern(phi[i,t]*z[i,t])
    }
  }
  # Priors
  a ~ dunif(-0.2,0.2)
  b ~ dunif(-0.2,0.2)
}

CJS_taille <- jags.parallel(data = jags.data,
                            inits = inits,
                            parameters.to.save = parameters,
                            model.file = cjs_taille,
                            n.chains = 4,
                            n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

#save(CJS_taille, file = "R/object/CJS_taille.RData")

####################################################################################################################
####################################################################################################################
####################################################################################################################


##################### Model 1
# avec b qui pourra varier avec les traitements pour phi

# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a = runif(1,-0.2,0.2), b = runif(1,-100,100), c = runif(1,-0.2,0.2), d = runif(1,0,100), z = zi)
}

parameters = c("a","b","c","d")

cjs_taille <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      phi[i,t_bis] <- 0
      p[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(phi[i,t]) <- a*s[i,t]+b
      logit(p[i,t]) <- c*s[i,t]+d
      z[i,t] ~ dbern(p[i,t]*z[i,t-1])
      y[i,t] ~ dbern(phi[i,t]*z[i,t])
    }
  }
  # Priors
  a ~ dunif(-0.2,0.2)
  b ~ dunif(-100,100)
  c ~ dunif(-0.2,0.2)
  d ~ dunif (-100,100)
}

CJS_taille <- jags.parallel(data = jags.data,
                            inits = inits,
                            parameters.to.save = parameters,
                            model.file = cjs_taille,
                            n.chains = 4,
                            n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

#save(CJS_taille, file = "R/object/CJS_taille.RData")


################ je comprends rien au logit

x = runif(1000,0,1)
y = logit(x)

df = data.frame(x,y)
ggplot(df)+
  geom_point(aes(x,y))


x = runif(1000,-10,10)
y = invlogit(x)

df = data.frame(x,y)
ggplot(df)+
  geom_point(aes(x,y))

x = runif(1000,-10,10)
y = 1/(1+exp(-x))

df = data.frame(x,y)
ggplot(df)+
  geom_point(aes(x,y))
