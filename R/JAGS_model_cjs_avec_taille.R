################# data
library(arules)
s <- s2

tr <- BDD_a %>% filter(Year != "2022") %>% group_by(Tag_id, Treatment) %>% summarize() %>% ungroup()
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

t <- data.frame(t$Tag_id, apply(t[2:7],2, function(x) ifelse(is.na(x),0,1)), t$Treatment, t$Lake)
t <- na.omit(t)

#set.seed(123)
#t <- t[sample(1:nrow(t), 1000),]

s <- s %>% filter (Tag_id %in% t$t.Tag_id)
t <- t %>% arrange(t.Tag_id)
s <- s %>% arrange(Tag_id)

s_group <- data.frame(apply(s[2:7], 2, function(x) discretize(x,method = "fixed",breaks = c(0,100,125,150,175,200,300), labels= c(1,2,3,4,5,6))))
s_group [] <- apply(s_group, 2, as.numeric)

CH <- as.matrix(t[2:7])

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
                  size = as.matrix(s[2:7]*10^-2),
                  size_group = s_group,
                  Treatment = t$t.Treatment,
                  f = apply(CH, 1, get.first), 
                  nind = nrow(CH), 
                  noccas = ncol(CH),
                  ni = 100000,
                  zi = z.inits(CH))



# premier chose trouver les bons interval pour les paramètres
# Ensuite partir avec des inits différents pour chaques chaines
# 


##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       a2_p = runif(1,-7,0), #a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(1,-20,10),
       z = zi)
}

#inits <- function(){
#  list(a_p = 8, a_phi = 4,
#       a2_p = -4, #a2_phi = -3.3,
#       b_p = -2, b_phi = -2.5,
#       z = zi)
#}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi") # , "a2_phi")

cjs_taille1 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*10^-4*size[i,t]^2 + a_p*10^-2*size[i,t] + b_p
      logit(phi[i,t]) <- a_phi*10^-2*size[i,t] + b_phi # + a2_phi*10^-4*size[i,t]^2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,13)
  a_phi ~ dunif(1,100)
  a2_p ~ dunif(-10,0)
  #a2_phi ~ dunif(-10,0)
  b_phi ~ dunif(-20,10)
  b_p ~ dunif(-10,5)
}


CJS_taille1 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille1,
                             n.chains = 3,
                             n.iter = ni)
#print(CJS_taille)
#autocorr.plot(CJS_taille, ask = F)
#traceplot(CJS_taille, ask = F)

save(CJS_taille1, file = "R/object/CJS_taille1.RData")

##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       a2_p = runif(1,-7,0), a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(1,-20,10),
       z = zi)
}

#inits <- function(){
#  list(a_p = 8, a_phi = 4,
#       a2_p = -4, a2_phi = -3.3,
#       b_p = -2, b_phi = -2.5,
#       z = zi)
#}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi", "a2_phi")

cjs_taille3 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*10^-4*size[i,t]^2 + a_p*10^-2*size[i,t] + b_p
      logit(phi[i,t]) <- a_phi*10^-2*size[i,t] + b_phi + a2_phi*10^-4*size[i,t]^2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,13)
  a_phi ~ dunif(1,100)
  a2_p ~ dunif(-10,0)
  a2_phi ~ dunif(-10,0)
  b_phi ~ dunif(-20,10)
  b_p ~ dunif(-10,5)
}


CJS_taille3 <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille3,
                             n.chains = 3,
                             n.iter = ni)

save(CJS_taille3, file = "R/object/CJS_taille3.RData")

##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       a2_p = runif(1,-7,0), a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(4,-20,10),
       z = zi)
}

#inits <- function(){
#  list(a_p = 8, a_phi = 4,
#       a2_p = -4, a2_phi = -3.3,
#       b_p = -2, b_phi = rep(-2.5,4),
#       z = zi)
#}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi", "a2_phi")

cjs_taille_treatment <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*10^-4*size[i,t]^2 + a_p*10^-2*size[i,t] + b_p
      logit(phi[i,t]) <- a_phi*10^-2*size[i,t] + b_phi[Treatment[i]] + a2_phi*10^-4*size[i,t]^2
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,13)
  a_phi ~ dunif(1,100)
  a2_p ~ dunif(-10,0)
  a2_phi ~ dunif(-10,0)
  
  b_p ~ dunif(-10,5)
  for (tr in 1:4){
    b_phi[tr] ~ dunif(-20,10)
  }
}


CJS_taille_treatment <- jags.parallel(data = jags.data,
                             inits = inits,
                             parameters.to.save = parameters,
                             model.file = cjs_taille_treatment,
                             n.chains = 3,
                             n.iter = ni)

##################### Model phi(~taille)p(~taille)


# logit(phi) <- / ~ c*taille +d
# logit(p) <- / ~ a*taille + b

inits <- function(){
  list(a_p = runif(1,5,13), a_phi = runif(1,1,20),
       a2_p = runif(1,-7,0), # a2_phi = runif(1,-7,0),
       b_p = runif(1,-10,5), b_phi = runif(4,-20,10),
       z = zi)
}

#inits <- function(){
#  list(a_p = 8, a_phi = 4,
#       a2_p = -4, a2_phi = -3.3,
#       b_p = -2, b_phi = rep(-2.5,4),
#       z = zi)
#}

parameters = c("a2_p","a_p","a_phi","b_p","b_phi") #, "a2_phi")

cjs_taille_treatment2 <- function(){
  # Likelihood
  for (i in 1:nind){
    z[i,f[i]] <- 1
    for (t_bis in 1:f[i]){
      p[i,t_bis] <- 0
      phi[i,t_bis] <- 0
    }
    for (t in (f[i]+1):noccas){
      logit(p[i,t]) <- a2_p*10^-4*size[i,t]^2 + a_p*10^-2*size[i,t] + b_p
      logit(phi[i,t]) <- b_phi[Treatment[i]] + a2_phi[Treatment[i]]*10^-4*size[i,t]^2 #a_phi[Treatment[i]]*10^-2*size[i,t] + 
      z[i,t] ~ dbern(phi[i,t]*z[i,t-1])
      y[i,t] ~ dbern(p[i,t]*z[i,t])
    }
  }
  # Priors
  a_p ~ dunif(5,13)
  a2_p ~ dunif(-10,0)
  b_p ~ dunif(-10,5)
  for (tr in 1:4){
    b_phi[tr] ~ dunif(-20,10)
    #a2_phi[tr] ~ dunif(-10,0)
    a_phi[tr] ~ dunif(1,100)
  }
}


CJS_taille_treatment2 <- jags.parallel(data = jags.data,
                                      inits = inits,
                                      parameters.to.save = parameters,
                                      model.file = cjs_taille_treatment2,
                                      n.chains = 3,
                                      n.iter = ni)
